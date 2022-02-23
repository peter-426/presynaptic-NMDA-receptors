%
% Main_STP_run_cpp.m
%

% This Matlab script uses the "system" command to run the STP model.
% The model is implemented in the AN_STP.cpp file, which is compiled by 
% by the GNU C++ compiler g++, or gcc depending on your platform.
% 
% The output of the C++ executable (a.out) is saved to csv files.
% The Matlab script AN_STP_plots.m loads the csv data files and 
% generates plots.

clear;
%close all;
clc;
%------------------
disp(' ... running ...');
s_model=1;

%------------------
if s_model == 1
   sensor_model = 'Hill';
elseif s_model==2
   sensor_model = 'Markov6';
elseif s_model==3
   sensor_model = 'Markov7';
elseif s_model==4
   sensor_model = 'Allosteric';
else
   return
end

linux = 1;  % else Windows

% from cmd line:$ g++ -w Main.cpp ./lib/*.h -DHill -I lib/ 

if (linux==1) 
    cmd_line=sprintf('g++ -w Main.cpp lib/*.cpp -D%s -I lib/', sensor_model);
    prog_cpp='./a.out';  
else 
    return
    %cmd_line=sprintf('c:\MinGW\bin\gcc.exe -w Main.cpp -D%s -I lib/', sensor_model);
    %prog_cpp='a.exe';
end

fprintf(cmd_line);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% status=system(cmd_line, '-echo');   % -w no warnings
% 
% if (status > 0)
%     fprintf('\n   compile time error  \n\n');
%     return;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


choices=7;  %[1,2,5,7];     % 1=1Hz,  2=5Hz,  3=10Hz,  4=75 isi,  5=20Hz,  6=40Hz,  7=50Hz   +astro+spine

TRIALS=10;

ap5=1;  ryr=0;  astro=1;        % <<<<<<<<<<<<<<<<<=====================

cond = sprintf('%d %d %d', ap5, ryr, astro);

fprintf('\n\nTrials=%d, ap5=%d,  ryr=%d,  astro=%d\n', TRIALS, ap5, ryr, astro);

cmd1 = sprintf('%s %d  %0.2f %d %s', prog_cpp, 1000, 10, TRIALS, cond);
cmd2 = sprintf('%s %d  %0.2f %d %s', prog_cpp, 200,   2, TRIALS, cond);
cmd3 = sprintf('%s %d  %0.2f %d %s', prog_cpp, 100,   1, TRIALS, cond);
cmd4 = sprintf('%s %d  %0.2f %d %s', prog_cpp, 75, 0.75, TRIALS, cond);
cmd5 = sprintf('%s %d  %0.2f %d %s', prog_cpp, 50, 0.50, TRIALS, cond);
cmd6 = sprintf('%s %d  %0.2f %d %s', prog_cpp, 25, 0.25, TRIALS, cond);
cmd7 = sprintf('%s %d  %0.2f %d %s', prog_cpp, 10, 0.10, TRIALS, cond);
cmd8 = sprintf('%s %d  %0.2f %d %s', prog_cpp, 75, 0.15, TRIALS, cond);

for ch=choices  % ch in list of menu choices

    ex='10-spikes';  % ex=experiment 
              
    % Test case designed to replicate Figure 10 of the paper by McGuinness et al.,
    % "Presynaptic NMDARs in the hippocampus facilitate transmitter release at 
    % theta frequency", 2010.
    
    if (ch == 100) 
        disp('Test 1,5, and 20 Hz');
        system(cmd1, '-echo');
        AN_STP_plots(1000, ex, sensor_model)
    
        system(cmd2, '-echo')
        AN_STP_plots(200, ex, sensor_model)   
    
        system(cmd5, '-echo');
        AN_STP_plots(50, ex, sensor_model)
    end
    
    
    % Test each case one at a time.
    %
    if (ch == 1)  % 1 Hz
        disp('1000 isi, 1 Hz');
        status = system(cmd1, '-echo');
        Main_STP_plots(1000, ex, sensor_model)
    end
    
    if (ch == 2)  % 5 Hz
        tic
        disp('200 isi, 5 Hz');
        system(cmd2, '-echo');
        Main_STP_plots(200, ex, sensor_model)  % <<< MAN_ vs AN_ and Linux is case sensitive
        toc    
    end
    
    if (ch == 3)  % 10 Hz
        disp('100 isi, 10 Hz');
        system(cmd3, '-echo');
        Main_STP_plots(100, ex, sensor_model)
    end
    
    if (ch == 4)  % 13.33 Hz, 75 ms ISI for PPF experiments. 
        disp('75 isi, 13.33 Hz');
        system(cmd4, '-echo');
        Main_STP_plots(75, ex, sensor_model)
    end
    
    if (ch == 5)   % 20 Hz
        system(cmd5, '-echo');
        Main_STP_plots(50, ex, sensor_model)
    end
    
    if (ch == 6)   % 40 Hz
        disp('25 isi, 40 Hz');
        system( cmd6 , '-echo');
        Main_STP_plots(25, ex, sensor_model)
        if astro==1
            Main_STP_plots_astro(25, sensor_model)
        end
    end
        
    if (ch == 7)   % 50 Hz
        disp('10 isi, 100 Hz');
        system( cmd7 , '-echo');
        Main_STP_plots(10, ex, sensor_model)
        
        if astro==1
            Main_STP_plots_astro(10, sensor_model)   % isi=20 ms
        end
    end
    
    % Ryanodine receptor experiments (Paired Pulse Facilitation).
    %
    % Test case designed to replicate results reported in the paper by
    % Emptage, Reid, and Fine, 
    % "Calcium stores in hippocampal synaptic boutons mediate short-term plasticity, 
    % store-operated Ca2+ entry, and spontaneous transmitter release", 2001.
    %
    if (ch == 200)   
        disp('75 isi, 13.33 Hz');
        ex='ppf';
        system(cmd8, '-echo');
        Main_STP_plots(75, ex, sensor_model)
    end

end % for loop

% 3 x 3 set of plots for 1,5, and 20 Hz
Main_STP_bar(sensor_model)



