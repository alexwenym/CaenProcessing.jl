
# properties of digitizer setup
const volts_per_ADC = 2u"V" / 2^14
const time_per_sample = 2u"ns"
const time_per_clock_cycle = 8ns #125 mHz clock cycle
const R_load = 50u"Ω" # load resistor

#const n_samples = 1000
const header_size = 6 

const T_HEADER = UInt32
const T_DATA = Int16

"""
    calc_n_evts( fname )
since we know the number of samples + the dtypes, 
we can pre_calc the number of events

NB: could add a check for divisibility? = check on whether file data is somehow non-standard
"""
calc_n_evts( fname, n_samples ) = filesize( fname ) ÷ ( 
    header_size * sizeof(T_HEADER) + n_samples * sizeof(T_DATA) )

# ====================================================

"""
    get_n_event_samples( fname )
Get the number of samples for every event in the file.
Since all events have the same number of samples,
just need to read the header of the first event
"""
function get_n_event_samples( fname )
    open( fname, "r" ) do io
        header = [ read(io, UInt32) for i in 1:6 ]
        return ( header[1] - 24 ) ÷ 2
    end
end

"""
    read_next_trigger!( io, waveform, event_id  )

Read the next trigger data from the binary input stream `io`,updating the values of `waveform` and `event_id`. 
"""
function read_next_trigger!( io, waveform::AbstractArray{Int16},
    event_id, trigger_time, num_overflows )

    try
        header = [ read(io, UInt32) for i in 1:6 ]

        event_size = ( header[1] - 24 ) ÷ 2 
        # board_id = header[2]
        # pattern = header[3]
        # channel = header[4]
        event_counter = header[5]

        trigger_time_tag = header[6]
        # problem: the digitizer saves "trigger_time_tag" as UInt32, which will overflow w/o warning.
        # to solve this, check if the time since prev. trigger is of the same order of mag. as typemax( Int32 ) ≃ 2e9
        t = convert(Float32, trigger_time_tag) + num_overflows[1] * typemax(Int32)
        # and then increment
        if ( trigger_time[1] - t ) > 2e9
            num_overflows[1] += 1
            t += typemax(Int32)
        end

        # if length(waveform) != event_size
        #     error("waveform vector has length $(length(waveform)), but the event length is $(event_size)!")
        # end

        read!(
            io,  
            waveform
        )
        event_id[1] = event_counter
        trigger_time[1] = t
        return true 

    catch error
        if ! isa( error, EOFError ) 
            display(error)
        end
        return false
    end
end

"""
    instantiate_data_vecs( n_samples )

Instantiate vectors for waveforms, event_id, trigger_time, etc. 
to be modified in `read_next_trigger!`. 
"""
function instantiate_data_vecs( n_samples )

    w = Vector{T_DATA}(undef, n_samples)
    event_id = Vector{T_HEADER}( undef, 1 ) 
    trigger_time = Vector{float(T_HEADER)}( undef, 1 )
    num_overflows = [0]

    return w, event_id, trigger_time, num_overflows
end

"""
    read_waveforms( f_dat; n_evts=Inf )
"""
function read_waveforms( f_dat, n_evts=Inf )

    n_samples = get_n_event_samples( f_dat )

    n_evts = ( isinf(n_evts) ) ? calc_n_evts(f_dat, n_samples) : n_evts 
    waveform_array = Matrix{T_DATA}(undef, n_samples, n_evts )

    # instantiate for while loop
    event_id, trigger_time, num_overflows = instantiate_data_vecs( n_samples )[2:end]

    open( f_dat, "r" ) do io     
    # read + process data evt-by-evt from file   
        for i in axes( waveform_array, 2 )
            read_next_trigger!( io, 
                (@view waveform_array[:, i]), 
                event_id, trigger_time, num_overflows )
        end
    end
    return waveform_array
end


begin
# calculation functions 

    function check_saturation( waveform )
        return any( waveform .== 0 )
    end

    function calc_baseline( waveform, xmin, xmax )
        mean( waveform[xmin:xmax] )
    end

    function calc_baseline_rms( waveform, xmin, xmax )
        rms = sqrt( mean( (waveform[xmin:xmax]).^2 ) )
    end


    """
        calc_CFD_sample( waveform, pct )

    Return the index of the last sample pre-peak for which the waveform is at pct of the peak amplitude. 
    Algorithm = constant fraction discriminator.
    """
    function calc_CFD_sample_index( waveform, pct )

        k = argmax( waveform )
        max_w = waveform[k]

        # go backwards from peak until amplitude = pct * max_w 

        for i in k:-1:1
            ( waveform[i] < max_w * pct ) && return i 
        end

        # no such sample found (?)
        return -1 
    end

end

# process a single waveform
function process_waveform( waveform; CFD_threshold=0.2 )

    # ...

        is_saturated = check_saturation( waveform )

        # calc baseline
        baseline_width = 59 # samples (1..60, stays well before pulse rising edge at ~sample 75)
        baseline_xmin = 1
        baseline_xmax = baseline_xmin + baseline_width
        baseline = calc_baseline( 
            waveform,
            baseline_xmin,
            baseline_xmax
        )

        # subtract baseline
        waveform .-= baseline
        waveform .*= -1     # flip negative-polarity (PMT) pulses to positive

        # calc baseline rms 
        baseline_rms = calc_baseline_rms( 
            waveform, 
            baseline_xmin, 
            baseline_xmax
        )

        
        # find max val
        k = argmax( waveform )
        max_t = (k-1)
        max_w = waveform[k]

        # compute pulse start time
        k_start = calc_CFD_sample_index( 
            waveform, CFD_threshold )
        event_time = (k_start - 1)

        # calculate charge integral of the pulse, window centered on peak
        int_pre  = 10   # samples before peak
        int_post = 40   # samples after peak
        x_min = max(1, k - int_pre)
        x_max = min(length(waveform), k + int_post)
        charge_integral = sum( waveform[x_min:x_max] )

        # calculate pre-pulse integral
        
        # compute pulse shape discrimination (PSD)
        # method 1. = choose short "fast" time window 

        # compute charge-weighted mean time of pulses
    
    return (
        is_saturated = is_saturated, 

        baseline = baseline * volts_per_ADC,
        baseline_rms = baseline_rms * volts_per_ADC,

        waveform_max = max_w * volts_per_ADC,
        waveform_max_time = max_t * time_per_sample,

        event_time_CFD = event_time * time_per_sample,
        charge_integral = charge_integral * volts_per_ADC / R_load * time_per_sample,
    )
end


"""
    process_data( f_dat; n_evts=Inf )

Read a .dat file named `f_dat` from the Caen digitizer to a Table object.
"""
function process_data( f_dat; n_evts=Inf )

    n_samples = get_n_event_samples( f_dat )
    n_evts = convert( Int, min( n_evts, calc_n_evts(f_dat, n_samples) ) )

    # instantiate for while loop
    waveform =  Vector{float(T_DATA)}(undef, n_samples)
    w, event_id, trigger_time, num_overflows = instantiate_data_vecs( n_samples )

    # properties to save 
    blank_col(T::Type) = Vector{T}(undef, n_evts)
    blank_col(T::Quantity) = StructArray{T}( undef, n_evts )
    blank_col(T::Type, u::FreeUnits ) = blank_col( typeof(one(T) * u) )
    t = Table(
        event_id = blank_col( Int32 ),
        trigger_time = blank_col( Float64, s ),
        is_saturated = blank_col( Bool ),
        baseline = blank_col( Float64, mV ),
        baseline_rms = blank_col( Float64, mV ),
        waveform_max_time = blank_col( Int32, ns ),
        waveform_max = blank_col( Float64, mV ),
        event_time_CFD = blank_col( Int32, ns ),
        charge_integral = blank_col( Float64, C),
    )

    open( f_dat, "r" ) do io     
    # read + process data evt-by-evt from file   
        for i in 1:n_evts

            not_done = read_next_trigger!( io, w, 
                event_id, trigger_time, num_overflows )

            # convert type for calculations
            waveform .= w 

            # calculations
            calcs = process_waveform( waveform )

            # copy to table
            begin
                t.event_id[i] = event_id[1]
                t.trigger_time[i] = trigger_time[1] * time_per_clock_cycle

                t.is_saturated[i] = calcs.is_saturated
                t.baseline[i] = calcs.baseline
                t.baseline_rms[i] = calcs.baseline_rms

                t.waveform_max[i] = calcs.waveform_max
                t.waveform_max_time[i] = calcs.waveform_max_time

                t.event_time_CFD[i]    = calcs.event_time_CFD
                t.charge_integral[i] = calcs.charge_integral
            end
        end
    end

    return t
end