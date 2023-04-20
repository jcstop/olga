%{

Example to plot .abf file.

abf = Axon Binary Format

Reader downloaded here:
https://se.mathworks.com/matlabcentral/fileexchange/22114-fcollman-abfload


2022-02-06 ejimsve
created
%}


function f()
    
    clc();
    
    fname = "C:\Users\olne0992\Desktop\21819010.abf";
    % fname = "c:\ejimsve\private\olga\21819010.abf";
    
    %%
    time_stop = 5*60;
    
    
    [d, si, h] = abfload(fname, "stop", time_stop);
    
    %% sampling frequency (Hz)
    % si = sampling interval in us
    fs = 1 / (si * 1e-6)
    
    %% extract used channel
    y = d(:, 1);
    
    %% number of samples
    N = numel(y);
    
    %% create time array
    t = (1/fs) * (0:N-1);
    
    %% plot
    tiledlayout(1, 1, "Tilespacing", "Compact", "Padding", "Compact");
    nexttile();
    
    plot(t, y);
    %plot(t, y, '.-');
    grid on;
    grid minor;
    
    xlabel("Time (sec)");
    ylabel("Voltage (mV)");
    
    xticks(gca(), 0 : 10 : time_stop);
    
    
    
end