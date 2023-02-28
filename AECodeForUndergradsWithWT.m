clear all;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This matlab code has been created for the analysis of Acoustic Emission data recorded in binary file type (.bin). To use this code:

%       1. Edit the source directory (source_dir) to refer to the folder in which your binary files are stored (remember to keep the backslash at the end).
%       2. Set the gain which was set on the signal conditioning unit (SCUGain) during the experiment. 
%       3. Set the gain which was set on the Preamplifier unit (PreAmpGain) during the experiment. 
%       4. Press "Run". 

% The code will plot graphs of the time-domain, frequency domain and Time-Frequency domain (Wavelet Transform) for each file.
% Please note that the wavelet transform plot is for a limited range of data in the time domain due to being computationally intensive. These can be saved manually from the plot window.
% Parameters of duration, rise-time, decay time, number of counts, peak ampitude, Peak Frequency, Frequency centroid, Weighted Peak Frequency and AE energy will be stored in the workspace. These can be plotted or exported as desired.
% This code can handle files with multiple AE channels, all plots will be labelled as to their corresponding file and channel.
% Variables in the workspace are arranged so each column corresponds to a file and each row corresponds to the channel.

% Code by Dr. A. Crawford - a.crawford2@rgu.ac.uk
% March 2021


%% Set Input parameters

source_dir='R:\R-PHD-AC-1103263\AE Data\3mm bar plb\Edge breaks\Point 1\'; % Change this to refer to the folder containing your .bin files. Ensure you retain the \ at the end of the line
SCUGain = 12;   % Gain at SCU
PreAmpGain= 20;     %gain at preamp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Files
source_files = dir(fullfile(source_dir, '*.bin'));
n=length(source_files);   %n is number of .bin files in each folder
for k=1:n; %Repeat for each data file
            fname=source_files(k).name;
            fname2=[source_dir fname];
            fid = fopen(fname2,'r','b');
            para=fread(fid,4,'float');
            fs=para(4);
            status=fseek(fid,16,'bof');
            a = fread(fid,inf,'float');
            fclose(fid);
            clear fid;
            ndata = length(a);
            x = (0:1/fs:(ndata-1)/(para(1)*fs))'; %Create time domain x-axis
         
 v=1;
for s=1:para(1):ndata
    q=1;
    for b=s:(s+para(1)-1)
        y(v,q)=a(b);
        q=q+1;
    end
    v=v+1;
end
clear s v q b;
for s = 1:para(1)	 
    ymean = mean(y(:,s));
    y(:,s) = y(:,s)-ymean;
    channels=y/((10^(PreAmpGain/20))*(10^(SCUGain/20)));   
end 
clear y ymean status a ans s
%% Analysis
    for channelno=1:para(1) %Repeat for each channel
    label = sprintf('%s%d%s%d', 'File: ',k,' Channel: ',channelno);

        
    % Time Domain Signal Parameters
    
    %Plot Time Domain Signal
    figure
    plot(x,channels(:,channelno),'k');
    ylabel('Amplitude (V)')
    xlabel('Time (sec)')
    title(label)
    grid on
    
    %Calculate Time Domain Parameters

    %Peak Amplitude
    [PeakAmp PeakLoc]=max(abs(channels(:,channelno)));
    PeakAmplitude(channelno,k)=PeakAmp;

    %Threshold (10% of peak amplitude)
    Threshold=0.1*PeakAmp;
    
    % Rise Time
    StartOfHit=find(abs(channels(:,channelno))>Threshold,1);
    RiseTime(channelno,k)=(PeakLoc-StartOfHit)/fs;
    
    % Decay Time
    EndOfHit=find(abs(channels(:,channelno))>Threshold,1,'last');
    DecayTime(channelno,k)=(EndOfHit-PeakLoc)/fs;
    
    % Duration
    Duration(channelno,k)=RiseTime(channelno,k)+DecayTime(channelno,k);
    
        %Energy
    Energy(channelno,k)=trapz(x(StartOfHit:EndOfHit),channels(StartOfHit:EndOfHit,channelno).^2);
    
    % No of counts
    noofcounts(channelno,k)=(length(find(diff((channels(StartOfHit:EndOfHit,channelno)-Threshold)>0)~=0)+1))/2;  

    % Frequency Domain Analysis
    %FFT
    [pyy1,f1] = pwelch(channels(StartOfHit:EndOfHit,channelno),[],[],[1000:1000:1000000],fs);
    FFT(k).FFT(channelno,:)=pyy1;
    figure
    plot(f1,pyy1,'k');
    ylabel('Power Spectral Density')
    xlabel('Frequency (Hz)')
    xlim([0 1000000])
    title(label)
    grid on
    
    %Peak Frequency
    [PeakFreqAmp PeakFreqLoc]=max(pyy1);
    PeakFreq(channelno,k)=f1(PeakFreqLoc);
    
    %Frequency Centroid
    FreqF=f1.*pyy1;
    IntFreqF=trapz(f1,FreqF);
    IntF=trapz(f1,pyy1);
    FreqCent(channelno,k)=IntFreqF/IntF;
    
    %Weighted Peak Frequency
    WeightedPeak(channelno,k)=sqrt(PeakFreq(channelno,k)*FreqCent(channelno,k));
    
    
    %Wavelet Transform
    WTPoints=2500;      %Set number of points to include in WT (Increase this too much and matlab will struggle)
    scale=2:0.5:110;
    wavelet=cwt(channels(1:WTPoints,channelno),scale,'morl');
   for n=1:length(scale);
    wavelet(n,:)=smooth(abs(wavelet(n,:)),10);  
    wavelet(n,:)=smooth(abs(wavelet(n,:)),10);  
   end
   F=scal2frq(scale,'morl',1/fs);
    
   %Create WT Plot
   figure
   contourf(x(1:WTPoints),F,abs(wavelet),100,'LineStyle','none')
   xlabel('Time (Sec)');
   ylabel('Frequency (Hz)');
   title(label) 
   colormap('jet');
   
%%   
%    Create combination plot 
scrsz = get(0,'ScreenSize');    %Get screen size
figure('Position',[scrsz(4)/2 scrsz(4)/10 scrsz(3)/2 scrsz(4)/2]);    %Create new figure based on screen size


        %Create PSD subplot
        ax1=subplot(2,4,1);
        plot(pyy1,f1,'k');
        ylim([min(F) max(F)]);
        ylabel('Frequency (Hz)');
        xlabel('PSD')
        grid on

        cx1=subplot(2,4,4);
        set(cx1,'Visible','off');

        %Create WT subplot
        ax2=subplot(2,4,2:3);
           contourf(x(1:WTPoints),F,abs(wavelet),100,'LineStyle','none')
            xlabel('Time (Sec)');
            ylabel('Frequency (Hz)');
        colormap('jet');
        title(label)
        %Add Colourbar
        c1=colorbar;
        c1.Position=[cx1.Position(1) cx1.Position(2) cx1.Position(3)/6 cx1.Position(4)];
        ylabel(c1,'WT Coeff.')

        xlim([0 x(WTPoints)]);
        ylim([min(F) max(F)]);

        %Add Time Domain Plot
        ax3=subplot(2,4,6:7);
            plot(x,channels(:,channelno),'k');
            ylabel('Amplitude (V)')
            xlabel('Time (sec)')
        xlim([0 x(WTPoints)]);
        ylabel('Amplitude (Norm.)');
        xlabel('Time (Sec)');
        grid on

        linkaxes([ax1,ax2],'y')
        linkaxes([ax2,ax3],'x')
        

%%
    end
end
clear channelno channels fname fname2 fs k n ndata para pyy1 SCUGain source_dir source_files PreAmpGain Threshold StartOfHit PeakLoc PeakAmp label EndOfHit x PeakFreqAmp PeakFreqLoc FreqF IntF IntFreqF WTPoints wavelet scrsz scale F f1

why