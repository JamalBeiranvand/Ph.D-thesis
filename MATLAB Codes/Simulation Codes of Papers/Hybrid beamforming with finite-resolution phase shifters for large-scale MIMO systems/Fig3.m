% All functions, that are used in this code, are described on the "MATLAb Functions" document.
 % you can see more detail in the functions codes. 
 
 % Author  : Jamal Beiranvand (Jamalbeiranvand@gmail.com)
 % google scholar: https://scholar.google.com/citations?user=S6LywwsAAAAJ&hl=en
 % Date : 14/05/2020
%%
clear;
clc
warning off
%% User pannel Parameters (you may have to change some values in the "Plot" section)
addpath('F:\Ph.D\Matlab\Functions') % Add path of the "Functions" folder.
b = 1;                      % Values of Bits, it can be a vector, you can use 'inf' for infinity resolution
Txz = 1;                      % Number of the transmitter antennas on z axis
Txy = 64;                      % Number of the transmitter antennas on y axis
Rxz = 1;                       % Number of the receiver antennas on z axis
Rxy = 16;                       % Number of the receiver antennas on y axis
Ns  = 4;                       % Number of data streams
NRF = 4;                      % Number of RF chains
Ncl =15;                       % Number of Channel clusters(Scatters)
Nray= 1;                       % Number of rays in each cluster
% AngSpread = 10/180*pi;       % Angle spread, standard deviation of the angles in azimuth and elevation both of Rx and Tx
SNR_dB = -10:2:10;             % Signal to noise ratio in dB
realization = 200;             % Iteration of simulation
%% Initialize Parameters
Tx  = [Txz,Txy];               % Srtucture of the transmit antanna array
Rx  = [Rxz,Rxy];               % Srtucture of the receiver antanna array 
Nt  = Txz*Txy;                 % Number of the transmit antennas
Nr  = Rxz*Rxy;                 % Number of the receive antennas
SNR = 10.^(SNR_dB./10);        
%% Allocate memory space for matrices
R_OMP=zeros(length(SNR_dB),realization);
R_opt=zeros(length(SNR_dB),realization);
R_shrb=zeros(length(SNR_dB),realization);
R_shrbq1=zeros(length(SNR_dB),realization);
R_shrbq2=zeros(length(SNR_dB),realization);
R_shrbq3=zeros(length(SNR_dB),realization);
%%
for r=1:length(SNR_dB)
    for reali=1:realization
        % generate channel
        [H,At,Ar,alpha]=Channel(Tx,Rx,Ncl);
        
        % Full Digital beamforming
        [W_Dig,F_Dig]=DigitalBeamforming(H,Ns);
        
        % Sohrabi Algorithm, Quantized after computing
        [WD,WRF,VRF,VD]=Sohrabi_Alg(H,NRF,Ns,SNR(r)); 
        Vt=VRF*VD;
        Wt=WRF*WD;
        
        VRF1=Qtz(VRF,2^1);
        WRF1=Qtz(WRF,2^1);
        Vt1=VRF1*VD;
        Wt1=WRF1*WD;
        
        VRF2=Qtz(VRF,2^2);
        WRF2=Qtz(WRF,2^2);
        Vt2=VRF2*VD;
        Wt2=WRF2*WD;
        
        VRF3=Qtz(VRF,2^3);
        WRF3=Qtz(WRF,2^3);
        Vt3=VRF3*VD;
        Wt3=WRF3*WD;

        % Sohrabi Algorithm, Quantized during computing
        [WDq1,WRFq1,VRFq1,VDq1]=Sohrabi_Alg(H,NRF,Ns,SNR(r),2^1); 
        Vtq1=VRFq1*VDq1;       
        Wtq1=WRFq1*WDq1;
       
        [WDq2,WRFq2,VRFq2,VDq2]=Sohrabi_Alg(H,NRF,Ns,SNR(r),2^2); 
        Vtq2=VRFq2*VDq2;       
        Wtq2=WRFq2*WDq2;
        
        [WDq3,WRFq3,VRFq3,VDq3]=Sohrabi_Alg(H,NRF,Ns,SNR(r),2^3); 
        Vtq3=VRFq3*VDq3;       
        Wtq3=WRFq3*WDq3;
        
        R_shrb(r,reali) = log2(det(eye(Nr) + Wt*inv((Wt)'*Wt)*(Wt)' * H * Vt * (Vt)' * H'));
        R_shrb1(r,reali) = log2(det(eye(Nr) + Wt1*inv((Wt1)'*Wt1)*(Wt1)' * H * Vt1 * (Vt1)' * H'));
        R_shrb2(r,reali) = log2(det(eye(Nr) + Wt2*inv((Wt2)'*Wt2)*(Wt2)' * H * Vt2 * (Vt2)' * H'));
        R_shrb3(r,reali) = log2(det(eye(Nr) + Wt3*inv((Wt3)'*Wt3)*(Wt3)' * H * Vt3 * (Vt3)' * H'));
        R_shrbq1(r,reali) = log2(det(eye(Nr) + Wtq1*inv((Wtq1)'*Wtq1)*(Wtq1)' * H * Vtq1 * (Vtq1)' * H'));
        R_shrbq2(r,reali) = log2(det(eye(Nr) + Wtq2*inv((Wtq2)'*Wtq2)*(Wtq2)' * H * Vtq2 * (Vtq2)' * H'));
        R_shrbq3(r,reali) = log2(det(eye(Nr) + Wtq3*inv((Wtq3)'*Wtq3)*(Wtq3)' * H * Vtq3 * (Vtq3)' * H'));

        R_opt(r,reali) = log2(det(eye(Ns) + SNR(r)/Ns * pinv(W_Dig) * H * F_Dig * (F_Dig)' * H' * W_Dig));

    end
end
R_shrb1(isnan(R_shrb1))=0;
R_shrb1(isinf(R_shrbq1))=0;

R_shrb2(isnan(R_shrb2))=0;
R_shrb2(isinf(R_shrb2))=0;

R_shrb3(isnan(R_shrb3))=0;
R_shrb3(isinf(R_shrb3))=0;

R_shrbq1(isnan(R_shrbq1))=0;
R_shrbq1(isinf(R_shrbq1))=0;

R_shrbq2(isnan(R_shrbq2))=0;
R_shrbq2(isinf(R_shrbq2))=0;

R_shrbq3(isnan(R_shrbq3))=0;
R_shrbq3(isinf(R_shrbq3))=0;

R_shrb(isnan(R_shrb))=0;
R_shrb(isinf(R_shrb))=0;

Xaxis=SNR_dB;
Data1= sum((R_opt),2)/realization;
Data2=sum(R_shrb,2)/realization;
Data3(:,1)=sum(R_shrbq1,2)./(realization-sum(R_shrbq1==0,2));
Data3(:,2)=sum(R_shrbq2,2)./(realization-sum(R_shrbq2==0,2));
Data3(:,3)=sum(R_shrbq3,2)./(realization-sum(R_shrbq3==0,2));
Data4(:,1)=sum(R_shrb1,2)./(realization-sum(R_shrb1==0,2));
Data4(:,2)=sum(R_shrb2,2)./(realization-sum(R_shrb2==0,2));
Data4(:,3)=sum(R_shrb3,2)./(realization-sum(R_shrb3==0,2));
% Data4= sum(R_OMP,2)/realization;
%% plot 
Fig=figure;
hold on
P1=plot(Xaxis,Data1);
P2=plot(Xaxis,Data2);
P3=plot(Xaxis,Data3);
P4=plot(Xaxis,Data4);
% P1 properties
set(P1,'LineStyle','-');           % '-' , '--' , ':' , '-.' , 'none'
set(P1,'Marker','none');       % 'o' , '+'  , '*' , 'x' ,'s' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram', 'none'
set(P1,'LineWidth',2);
set(P1,'Color','k');
set(P1,'MarkerSize',7);
P1_group= hggroup;set(P1,'Parent',P1_group);
set(get(get(P1_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on'); 
% P1.MarkerIndices =[ 1 6];

% P2 properties
set(P2,'LineStyle',':');           % '-' , '--' , ':' , '-.'
set(P2,'Marker','none');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
set(P2,'LineWidth',2);
set(P2,'Color','k');
set(P2,'MarkerSize',7);
P2_group= hggroup;set(P2,'Parent',P2_group);
set(get(get(P2_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');

% P3 properties
set(P3,'LineStyle','-');           % '-' , '--' , ':' , '-.'
set(P3,'Marker','o');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
set(P3,'LineWidth',2);
set(P3,'Color','b');
set(P3,'MarkerSize',7);
P3_group= hggroup;set(P3,'Parent',P3_group);
set(get(get(P3_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');

% P3 properties
set(P4,'LineStyle','-');           % '-' , '--' , ':' , '-.'
set(P4,'Marker','s');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
set(P4,'LineWidth',2);
set(P4,'Color','r');
set(P4,'MarkerSize',7);
P4_group= hggroup;set(P4,'Parent',P4_group);
set(get(get(P4_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');

lgd=legend('Optimal Fully Digital BF','Proposed Alg. for $b=\infty$',' Proposed Alg. for $b$','Quantized $-$ Proposed Alg. for $b=\infty$');

% Axes Properties 
ax=gca;
ax.XLim = [min(Xaxis) max(Xaxis)];
ax.YLim = [3 36];
ax.XTick = Xaxis;
ax.XTickLabel=Xaxis;
ax.XLabel.String='SNR (dB)';
ax.YLabel.String='Spectral Efficiency (bits/sec/Hz)';
ax.XGrid='on';
ax.YGrid='on';
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'on';
ax.TickLabelInterpreter='latex';
ax.YLabel.Interpreter='latex';
ax.XLabel.Interpreter='latex';
ax.FontSize = 11;
ax.LineWidth = .9;
ax.Box='on';

% Legend Properties %
lgd.Location= 'northwest'; %%
lgd.Box='on';
lgd.LineWidth=.9;
lgd.Interpreter='latex';
lgd.FontSize = 10;

% Figure Properties
ax.Units='centimeters';
ax.Position = [ax.TightInset(1) ax.TightInset(2) 11 8];
Fig.Units='centimeters';
Fig.Position = [10 5 (ax.TightInset(1)+ ax.TightInset(3)+ ax.Position(3) ) (ax.TightInset(2)+ ax.TightInset(4)+ ax.Position(4))];

Ell=annotation('ellipse');
X_points=[1 2];
Y_points=[14 18];
Ell.Units='centimeters';
Dx=ax.Position(3)/(diff(ax.XLim));
Dy=ax.Position(4)/(diff(ax.YLim));
x_begin=ax.Position(1)-ax.XLim(1)*Dx+X_points(1)*Dx;
y_begin=ax.Position(2)-ax.YLim(1)*Dy+Y_points(1)*Dy;
x_End=(ax.Position(1)-ax.XLim(1)*Dx+X_points(2)*Dx)-x_begin;
y_End=(ax.Position(2)-ax.YLim(1)*Dy+Y_points(2)*Dy)-y_begin;
Ell.Position=[x_begin y_begin  x_End y_End];
% 
Textarrow1=annotation('textarrow');
X_points=[3 2];
Y_points=[13 16];
str1='$b=1$';
Textarrow1.String={str1};
Textarrow1.Interpreter='latex';
Textarrow1.FontSize = 10;
Textarrow1.Units='centimeters';
Dx=ax.Position(3)/(diff(ax.XLim));
Dy=ax.Position(4)/(diff(ax.YLim));
x_begin=ax.Position(1)-ax.XLim(1)*Dx+X_points(1)*Dx;
y_begin=ax.Position(2)-ax.YLim(1)*Dy+Y_points(1)*Dy;
x_End=(ax.Position(1)-ax.XLim(1)*Dx+X_points(2)*Dx)-x_begin;
y_End=(ax.Position(2)-ax.YLim(1)*Dy+Y_points(2)*Dy)-y_begin;
Textarrow1.Position=[x_begin y_begin  x_End y_End];

% %%
Ell2=annotation('ellipse');
X_points=[2 3];
Y_points=[21 23];
Ell2.Units='centimeters';
Dx=ax.Position(3)/(diff(ax.XLim));
Dy=ax.Position(4)/(diff(ax.YLim));
x_begin=ax.Position(1)-ax.XLim(1)*Dx+X_points(1)*Dx;
y_begin=ax.Position(2)-ax.YLim(1)*Dy+Y_points(1)*Dy;
x_End=(ax.Position(1)-ax.XLim(1)*Dx+X_points(2)*Dx)-x_begin;
y_End=(ax.Position(2)-ax.YLim(1)*Dy+Y_points(2)*Dy)-y_begin;
Ell2.Position=[x_begin y_begin  x_End y_End];

Textarrow2=annotation('textarrow');
X_points=[5 3];
Y_points=[15 22];
str1='$b=2$';
Textarrow2.String={str1};
Textarrow2.Interpreter='latex';
Textarrow2.FontSize = 10;
Textarrow2.Units='centimeters';
Dx=ax.Position(3)/(diff(ax.XLim));
Dy=ax.Position(4)/(diff(ax.YLim));
x_begin=ax.Position(1)-ax.XLim(1)*Dx+X_points(1)*Dx;
y_begin=ax.Position(2)-ax.YLim(1)*Dy+Y_points(1)*Dy;
x_End=(ax.Position(1)-ax.XLim(1)*Dx+X_points(2)*Dx)-x_begin;
y_End=(ax.Position(2)-ax.YLim(1)*Dy+Y_points(2)*Dy)-y_begin;
Textarrow2.Position=[x_begin y_begin  x_End y_End];

Ell3=annotation('ellipse');
X_points=[5 5.5];
Y_points=[27 28];
Ell3.Units='centimeters';
Dx=ax.Position(3)/(diff(ax.XLim));
Dy=ax.Position(4)/(diff(ax.YLim));
x_begin=ax.Position(1)-ax.XLim(1)*Dx+X_points(1)*Dx;
y_begin=ax.Position(2)-ax.YLim(1)*Dy+Y_points(1)*Dy;
x_End=(ax.Position(1)-ax.XLim(1)*Dx+X_points(2)*Dx)-x_begin;
y_End=(ax.Position(2)-ax.YLim(1)*Dy+Y_points(2)*Dy)-y_begin;
Ell3.Position=[x_begin y_begin  x_End y_End];

Textarrow3=annotation('textarrow');
X_points=[7 5.5];
Y_points=[17 27];
str1='$b=3$';
Textarrow3.String={str1};
Textarrow3.Interpreter='latex';
Textarrow3.FontSize = 10;
Textarrow3.Units='centimeters';
Dx=ax.Position(3)/(diff(ax.XLim));
Dy=ax.Position(4)/(diff(ax.YLim));
x_begin=ax.Position(1)-ax.XLim(1)*Dx+X_points(1)*Dx;
y_begin=ax.Position(2)-ax.YLim(1)*Dy+Y_points(1)*Dy;
x_End=(ax.Position(1)-ax.XLim(1)*Dx+X_points(2)*Dx)-x_begin;
y_End=(ax.Position(2)-ax.YLim(1)*Dy+Y_points(2)*Dy)-y_begin;
Textarrow3.Position=[x_begin y_begin  x_End y_End];

