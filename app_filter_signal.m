%{

Processing steps
 1. Decimate signal (reduce sample rate).
 2. Filter signal.


2022-02-06 ejimsve
created
%}


function app_filter_signal()
    
    clc();
    
    % desired low-pass cutoff frequency
    p.lowpass_fc = 0.5;
    
    % file analyze
    %p.fname = "c:\ejimsve\private\olga\21819010.abf";
    %p.fname = "C:\Users\olne0992\Desktop\21819010.abf";
    p.fname = "d:\olga\abf\21n09002.abf";
    
        
    % analyze segment 1
    p.time_start  =  0;
    p.time_length = 35;
    p.segment_nr  = 1;
    p.channel_nr  = 3;
    s1 = analyze_segment(p);
    
    % analyze segment 2
    p.time_start  = 460;
    p.time_length = 90;
    p.segment_nr  = 2;
    p.channel_nr  = 3;
    s2 = analyze_segment(p);
    
    
    
    %% compare filtered signal
    figure(900);
    tiledlayout(2, 1, "Tilespacing", "Compact", "Padding", "Compact");

    nexttile();
    plot(s1.t_filt, s1.x_filt);
    grid on;
    title("Segment 1");
    xlabel("Time (sec)");
    
    nexttile();
    plot(s2.t_filt, s2.x_filt);
    grid on;
    title("Segment 2");

    set(gcf(), "Name", "Compare filtered");

    
    %% compare derivative signal
    figure(901);
    tiledlayout(2, 1, "Tilespacing", "Compact", "Padding", "Compact");

    nexttile();
    plot(s1.t_diff, s1.x_diff);
    grid on;
    title("Segment 1");
    xlabel("Time (sec)");
    ylabel("Derivative (mV/s)")
    
    nexttile();
    plot(s2.t_diff, s2.x_diff);
    grid on;
    title("Segment 2");

    set(gcf(), "Name", "Compare derivative");
    
    
    
    
    %return
    
    %% Compare derivative max values
    figure(800);
    tiledlayout(1, 1, "Tilespacing", "Compact", "Padding", "Compact");
    nexttile();
    
    hold off;
    plot(s1.diff_max_idxes, s1.diff_max_values, '.-', "color", [0 0 1], "DisplayName", "Segment 1");
    
    hold on;
    plot(s2.diff_max_idxes, s2.diff_max_values, '.-', "color", [1 0 0], "DisplayName", "Segment 2");
    
    grid on;
    title("Derivative max values");
    xlabel("Peak index");
    ylabel("Max dv/dt (mV/s)");
    legend();

    set(gcf(), "Name", "Derivative max comparison");
    
    %% compare derivative min values
    figure(810);
    tiledlayout(1, 1, "Tilespacing", "Compact", "Padding", "Compact");
    nexttile();
    
    hold off;
    plot(s1.diff_min_idxes, s1.diff_min_values, '.-', "color", [0 0 1], "DisplayName", "Segment 1");
    hold on;
    plot(s1.diff_min_idxes([1 end]), s1.diff_min_avg * [1 1], '--', "color", [0 0 1], "DisplayName", "Segment 1 (avg)");
    
    
    hold on;
    plot(s2.diff_min_idxes, s2.diff_min_values, '.-', "color", [1 0 0], "DisplayName", "Segment 2");
    plot(s2.diff_min_idxes([1 end]), s2.diff_min_avg * [1 1], '--', "color", [1 0 0], "DisplayName", "Segment 2 (avg)");
    s2
    
    grid on;
    title("Derivative min values");
    xlabel("Peak index");
    ylabel("Max dv/dt (mV/s)");
    legend();

    set(gcf(), "Name", "Derivative min comparison");
    
    fprintf("Done!\n")
end



function o = analyze_segment(p)
    
    time_start  = p.time_start;
    time_length = p.time_length;
        
    
    %% 
    time_stop = time_start + time_length;
    
    % load file
    [d, si, h] = abfload(p.fname, "start", time_start, "stop", time_stop);
    fprintf("\n");
    drawnow();
    
    %% sampling frequency (Hz)
    % si = sampling interval in us
    fs = 1 / (si * 1e-6);
    
    %% extract used channel
    x_raw = d(:, p.channel_nr);
    
    %% number of samples
    N_raw = numel(x_raw);
    
    %% create time array
    t_raw = (1/fs) * (0:N_raw-1).';
    
    
    %% generate known signal for testing
    if 0
        A = 10;
        f = 0.3;
        x_raw = -40 + A * sin(2*pi*f*t_raw);
        
        x_prim_max_theory = A * 2*pi*f
    end
    
    % keep track of total delay for decimated and filtered signal
    delay = 0;
    
    % x_dec will be the decimated signal
    x_dec = x_raw;
    
    % fs2 will be the sample rate after decimation
    fs2 = fs;
    for k = 1 : 4
        
        % decimation factor
        R = 5;
        
        x_dec = decimate(x_dec, R, "FIR");
        %delay = delay  + 1/fs*30/2;
        fs2 = fs2/R;
    end
    
    %diff_dec_max = max(diff(x_dec) * fs2)
    
    %% time for decimated samples
    N_dec = numel(x_dec);
    t_dec = (1/fs2) * (0:N_dec-1).';
    t_dec = t_dec - delay;
    
        
    %% do final low-pass filtering
    if 1
        % order of FIR filter
        M = 100;
        
        % design FIR filter for filtering
        b = fir1(M, 2*p.lowpass_fc/fs2);
        
        % do the filtering
        x_filt = filter(b, 1, x_dec);
        
        % keep track of total delay
        delay = delay + (M/2) * 1/fs2;
    end
    
    %diff_filt_max = max(diff(x_filt(100:end)) * fs2)
    
    %% time for filtered samples
    N_filt = numel(x_filt);
    t_filt = (1/fs2) * (0:N_filt-1).';
    t_filt = t_filt - delay;
    
    
    %% differentiator to calc deriviative
    Nf    = 50; 
    Fpass = 0.1;
    Fstop = 4;

    % calculate coefficients for derivative filter
    d = designfilt('differentiatorfir', ...
        'FilterOrder', Nf, ...
        'PassbandFrequency', Fpass, ...
        'StopbandFrequency', Fstop, ...
        'SampleRate', fs2);
    
    % do the filtering to estimate derivative
    x_diff = filter(d, x_filt) * (fs2);
    
    % calculate delay for derivative signal
    delay_diff = d.FilterOrder / 2 * (1/fs2);
    
    % time difference between filtered signal and derivative
    t_diff = t_filt - delay_diff;
    
    
    %% remove transient
    
    % number of initial samples to skip
    num_trans = 100;
    
    % remove the samples
    t_diff = t_diff(1 + num_trans : end);
    x_diff = x_diff(1 + num_trans : end);
    
    
    %diff_diff_max = max(x_diff(100:end))
    
    %% print some info
    fprintf("Original sample rate = %.0f Hz\n", fs);
    fprintf("Decimated sample rate = %.0f Hz\n", fs2);
    


    %%
    %[filt_max_values, filt_max_times] = findpeaks(+x_filt, t_filt, "MinPeakProminence", 10);
    %x_filt_at_max = interp1(t_filt, x_filt, diff_max_times);



    
    %% find peak of derivative
    [diff_max_values, diff_max_times] = findpeaks(+x_diff, t_diff, "MinPeakProminence", 10);
    [diff_min_values, diff_min_times] = findpeaks(-x_diff, t_diff, "MinPeakProminence", 10);
    diff_min_values = -diff_min_values;
    
    %
    x_dec_at_max = interp1(t_filt, x_filt, diff_max_times);
    x_dec_at_min = interp1(t_filt, x_filt, diff_min_times);
    
    %
    num_max = numel(diff_max_values);
    diff_max_idxes = (1:num_max).';
    diff_max_avg = mean(diff_max_values);
    
    %
    num_min = numel(diff_min_values);
    diff_min_idxes = (1:num_min).';
    diff_min_avg = mean(diff_min_values);
    
    %% outputs
    o.diff_max_idxes  = diff_max_idxes;
    o.diff_max_values = diff_max_values;
    o.diff_max_avg    = diff_max_avg;
    
    o.diff_min_idxes  = diff_min_idxes;
    o.diff_min_values = diff_min_values;
    o.diff_min_avg    = diff_min_avg;
    
    o.t_filt = t_filt;
    o.x_filt = x_filt;
    
    o.t_diff = t_diff;
    o.x_diff = x_diff;
    
    
    %% plot raw + decimated
    figure(100 + p.segment_nr);
    tiledlayout(1, 1, "Tilespacing", "Compact", "Padding", "Compact");
    nexttile();
    
    % raw data
    hold off;
    plot(t_raw, x_raw, "DisplayName", "Original Data");
    
    % decimated data
    hold on;
    plot(t_dec, x_dec, "color", [0 1 1], "DisplayName", "Decimated");
    
    % filtered data
    hold on
    plot(t_filt, x_filt, "color", [1 0 0], "linewidth", 1, "DisplayName", "Filtered");
    
    % peak derivative
    hold on;
    plot(diff_max_times, x_dec_at_max, '.', "color", [1 0 0], "markersize", 20, "DisplayName", "Max slope");
    plot(diff_min_times, x_dec_at_min, '.', "color", [0 1 0], "markersize", 20, "DisplayName", "Min slope");
        
    %
    grid on;
    grid minor;
    xlim([0 time_length]);
    
    xlabel("Time (sec)");
    ylabel("Voltage (mV)");
    legend();
    title(sprintf("Raw and decimated data (segment %d)", p.segment_nr));
    f = gcf(); 
    f.Name = "Raw and decimated";

    
    
    %% plot decimated + filtered
    figure(200 + p.segment_nr);
    tiledlayout(1, 1, "Tilespacing", "Compact", "Padding", "Compact");
    nexttile();
    
    % decimated
    hold off;
    %plot(t_dec, x_dec);
    
    hold off
    plot(t_filt, x_filt, "linewidth", 1, "DisplayName", "Filtered");
    ylim([-80 -30]);
    %ylim([-100 0]);
            
    grid on;
    grid minor;
    xlim([0 time_length]);
    
    xlabel("Time (sec)");
    ylabel("Voltage (mV)");
    legend();
    
    str = sprintf("Decimated and filtered data, fc = %.2f Hz (segment %d)", p.lowpass_fc, p.segment_nr);
    title(str);
    set(gcf(), "Name", "Decimated and filtered");

    %% plot peaks
    %hold on
    %plot(filt_max_times, x_dec_at_max, '.', "color", [1 0 0], "markersize", 20, "DisplayName", "Max slope");


        
    
    %% plot derivative
    figure(300);
    tiledlayout(1, 1, "Tilespacing", "Compact", "Padding", "Compact");
    nexttile();
    
    hold off;
    plot(t_diff, x_diff);
    ylim(20 * [-1 1]);
        
    grid on;
    grid minor;
    xlim([0 time_length]);
    
    xlabel("Time (sec)");
    ylabel("Voltage derivative (mV/s)");
    title("Derivative of filtered signal");

    set(gcf(), "Name", "Derivate of filtered");

    
    
    
    %% plot peak derivative values
    figure(500);
    num_max = numel(diff_max_values);
    diff_max_idxes = (1:num_max).';
    
    diff_max_avg = mean(diff_max_values);
    
    hold off;
    plot(diff_max_idxes, diff_max_values, '.-', "markersize", 10);
    
    hold on;
    plot(diff_max_idxes([1 end]), diff_max_avg * [1 1], '.-');
    
    xlabel("Peak index");
    ylabel("Max dv/dt (mV/s)");
    grid on;
    grid minor;
    set(gcf(), "Name", "Peak derivative values");


    %% write csv file
    data = [diff_max_idxes diff_max_values];
    fname = sprintf("derivative_max_segment_%d.csv", p.segment_nr);
    writematrix(data, fname);
    
    
    %% plot dip derivative values
    figure(510);
    num_min = numel(diff_min_values);
    diff_min_idxes = (1:num_min).';
    
    diff_min_avg = mean(diff_min_values);
    
    hold off;
    plot(diff_min_idxes, diff_min_values, '.-', "markersize", 10);
    
    hold on;
    plot(diff_min_idxes([1 end]), diff_min_avg * [1 1], '.-');
    
    xlabel("Peak index");
    ylabel("Min dv/dt (mV/s)");
    grid on;
    grid minor;
    set(gcf(), "Name", "Dip derivative values");
    
    
    %% write csv file
    data = [diff_min_idxes diff_min_values];
    fname = sprintf("derivative_min_segment_%d.csv", p.segment_nr);
    writematrix(data, fname);
   
  
    
    
    
    
end