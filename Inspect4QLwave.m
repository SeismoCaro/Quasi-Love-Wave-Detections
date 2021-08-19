%% Inspect4QLwave written by Caroline Eakin ANU January 2020

% Matlab function to plot surface waves and look for Quasi-Love Waves
% Requires visual inspection of waveforms
% Will ask to either keep (1) or discard (0) the event 

% Code reproduces Figure S1 of following paper:
% Eakin (202-), Communications Earth & Environment
% Pre-print: https://www.researchsquare.com/article/rs-121788/v1 

%% %%%% INPUTS %%%%%
% Load example file "TestData_MORW_FigS1.mat" which contains event info
% from Figure S1

% E = East component 
% N = North component
% Z = Vertical component
% t = time
% sdt = sampling rate 
% dis = epicentral distance of event
% bazi = event backazimuth
% slat = station latitude
% slon = station longitude 

%% %%%% OUTPUTS %%%%%
% column 1 = delay time between QL and max G1 (dt) i.e. cross-correlation lag time
% column 2 = fractional distance along ray path of scatterer (dx)
% column 3 = QL scatterer latitude
% column 4 = QL scattter longitude
% column 5 = value of cross-correlation between QL and G1 (could be either + or -)
% column 6 = quality assessment by user, either "1" to keep or "0" discard
% column 7 = SNR of G1 compared to beginning of seismogram

%%

function out = Inspect4QLwave(E,N,Z,t,sdt,dis,bazi,slat,slong)

spos=get(0,'ScreenSize');
width= spos(3)-10 ; height=spos(4)-10;
xpos = spos(1) +5;  ypos = height - 100 ;
figpos =[xpos ypos width height];

si=1; % 1 event

dtQl=[]; dx=[]; scatlat=[]; scatlon=[]; qualtable=[];
 

% filter
f1 = 0; 
f2 = 1/100; % Hz 

ny = round(1 / sdt / 2); %nyquist frequency
[b,a]  = butter(3, [f2]/ny,'low');
Zf = (filtfilt(b, a, Z));
Ef = (filtfilt(b, a, E));
Nf = (filtfilt(b, a, N));

% Rotate to LQT
inc=0;
bazi2 = bazi/180*pi; 

M = [cos(inc)     -sin(inc)*sin(bazi2)    -sin(inc)*cos(bazi2);
     sin(inc)      cos(inc)*sin(bazi2)     cos(inc)*cos(bazi2);
        0              -cos(bazi2)             sin(bazi2)];

% M = rot3D(inc, bazi);
ZEN = [Zf, Ef, Nf]';         %
LQT = M * ZEN; %rotating    

Lo = LQT(1,:);                                       %
Qo = LQT(2,:);                                       %
To = LQT(3,:);  

scal=1;
L=scal*Lo/max(abs(Lo));
Q=scal*Qo/max(abs(Qo));
T=scal*To/max(abs(To));

% plot relative to max Love arrival
[G1max pos]=max(abs(T));
maxt=pos*sdt;
[R1max posR]=max(abs(L));

% calculate SNR 
Tnoise=max(abs(T(1:(300/sdt))));
Tsig=G1max;
snr=Tsig/Tnoise;

if snr>5
%     si=si+1; % Use if want to scale up to more than 1 event
    
    % Plot each seismogram
    xj=figure('Position',figpos);
    set(gcf,'color','w')
% %     eqtitle=strcat('EQ: ',num2str(j),'  /    ', eq(j).dstr,' / SNR: ',num2str(snr),' / Mw: ',num2str(roundn(eq(j).Mw,-1)), ' / Dis:  ',num2str(round(eq(j).dis)), ' /  BA: ', num2str(round(eq(j).bazi)),' / ',eq(j).region);
    
    subplot(611)
    plot(t,T,'color',[0 0.5 0],'LineWidth',1.5)
    ylabel('T; Transverse')
    hold all
    plot([t(pos) t(pos)],[-1 1],'--','color',[0 0.5 0],'LineWidth',1.0)
    grid on
% %     title(eqtitle)
    
    subplot(612)
    plot(t,L,'b','LineWidth',1.5); hold all
    plot([t(pos) t(pos)],[-1 1],'--','color',[0 0.5 0],'LineWidth',1.0)
    plot([t(posR) t(posR)],[-1 1],'--','color','b','LineWidth',1.0)
    ylabel('L: Vertical')
    hold all
    grid on
    
    subplot(613)
    plot(t,Q,'r','LineWidth',1.5); hold all
    plot([t(pos) t(pos)],[-1 1],'--','color',[0 0.5 0],'LineWidth',1.0)
    ylabel('Q: Radial')
    grid on
    
    
    subplot(614)
    plot(t,L,'b','LineWidth',1.5)
    hold all
    Qgrad=gradient(Q)./max(abs(gradient(Q)));
    plot(t,Qgrad,'r--','LineWidth',1.5);
    grid on
    plot([t(pos) t(pos)],[-1 1],'--','color',[0 0.5 0],'LineWidth',1.0)
    ylabel('L & gradient Q')
    axis([0 5000 -1 1])
    
    offs=1000/sdt;
    Q3=zeros(1,length(L)+(2*offs));
    Q3(offs:offs-1+length(Q))=Q;
 
    
    subplot(615)
    LQg=L+Qgrad;
    plot(t,LQg,'m','LineWidth',1.5); hold all
    plot([t(pos) t(pos)],[-2 2],'--','color',[0 0.5 0],'LineWidth',1.0)
    grid on
    ylabel('L+Qgrad stack')


%% Cross-correlation (Xcorr) to find QL wave

[r, lags]=xcorr(LQg,T);

ind=[]; ii=0;
gr=gradient(abs(r)); 
    for jr=1:length(gr)-1
        if (gr(jr)>0)&&(gr(jr+1)<0)
            ii=ii+1; ind(ii)=jr;
        end
    end
    
    ii=0; ind2=0;
    for jr=2:length(ind)-1
        if (abs(r(ind(jr)))>abs(r(ind(jr-1))))&&(abs(r(ind(jr)))>abs(r(ind(jr+1))))
            ii=ii+1; ind2(ii)=jr;
        end
    end
    
    qli=find(lags(ind(ind2))>0,1,'first');
    qldp=lags(ind(ind2(qli)));
    qldpt=(qldp*sdt);
    
    if isempty(qldpt)==1
        close(gcf)
    end
    
    subplot(615)
    if r(ind(ind2(qli)))>0 % positive correlation
        plot(t+qldpt,T,'color',[0 0.5 0],'LineWidth',1.5)

    else % negative correlation
        plot(t+qldpt,-T,'color',[0 0.5 0],'LineWidth',1.5)
    end
    try
        plot([t(pos+qldp) t(pos+qldp)],[- 2 2],':','color',[0 0.5 0],'LineWidth',1.5)
    catch
        close(gcf)
    end
    xlim([0 5000])
    

subplot(616)
plot(sdt*lags,abs(r),'k'); hold all; 
plot(sdt*lags(ind),abs(r(ind)),'x-'); 
plot(sdt*lags(ind(ind2)),abs(r(ind(ind2))),'bo')
plot(sdt*lags(ind(ind2(qli))),abs(r(ind(ind2(qli)))),'g*')
plot([sdt*lags(ind(ind2(qli))) sdt*lags(ind(ind2(qli)))],[abs(r(ind(ind2(qli)))) 5000],':','color',[0 0.5 0],'LineWidth',1.5)
axis([0-maxt 5000-maxt 0 5000])
ylabel('Xcorr')
plot([0 0],[0 5000],'--','color',[0 0.5 0],'LineWidth',1.0)
grid on
        

%% User Input & calculate distance
pause 
close(gcf)
userqual = input('Keep (1), or Discard (0)? : ')

dtQL(si)=sdt*lags(ind(ind2(qli)));

G1time=t(pos);
R1time=t(posR);

dx(si)=(dtQL(si)*dis)/(R1time-G1time);

[scatlat(si) scatlon(si)] = reckon(slat,slong,dx(si),bazi);

out(si,1:7)=[dtQL(si) dx(si) scatlat(si) scatlon(si) r(ind(ind2(qli))) userqual snr];

end

end










