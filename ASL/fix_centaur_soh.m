function w = fix_centaur_soh(w, startdate)
    SECS_PER_DAY = 60*60*24;
    for c=1:numel(w)
        chan = get(w(c),'channel');
        samprate = get(w(c),'freq');
        nsamps = get(w(c),'data_length');
        snum = get(w(c),'start');
        enum = get(w(c),'end');
        if nsamps == 0
            enum = snum;
            samprate = NaN;
        else
            if snum>0 & enum>0
                duration_in_seconds = (enum-snum)*SECS_PER_DAY;
            else
                duration_in_seconds = SECS_PER_DAY;
            end
            if samprate == 0
                samprate = nsamps/duration_in_seconds;
            end        
            if snum>0 & enum==0
                enum = (snum + 1) - 1/(samprate*SECS_PER_DAY);
            end
        end
        w(c) = set(w(c),'freq',samprate);    
            
end