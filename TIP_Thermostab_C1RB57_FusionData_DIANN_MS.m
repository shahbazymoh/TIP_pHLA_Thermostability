%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% A pipeline for thermostability profiling of HLA-bound peptides %%%%%%
%%%%%%%%% Thermal immunopeptidome profiler (TIP) code, version 1.2 %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% This code has been created, developed, and validated by Mohammad (Moh) Shahbazy 
%%% A PhD canidate at Department of Biochemistry and Molecular Biology  
%%% Clayton Campus, Monash University (Melbourne, Australia)
%%% Laboratory: Purcell Lab, Biomedicine Discovery Institute (BDI) 
%%% Acknowledgements: Dr Nathan Croft and Prof Anthony Purcell, for supervision
%%% Copyright© Mohammad Shahbazy (2019-2023).



%%% Creation date: June 2019
%%% Last modification: March 2023

% EXAMPLE DATASET: C1R-B*57:01 data acquired by Fusion Orbitrap (Thermo)

% DIA-NN version 1.8.0 was used for library-based DIA data processing
% Compatible with the next versions of DIA-NN (e.g., 1.8.1 or 1.9.2)

%%%%%%%%%%%%%%%*************** IMPORTANT ***************%%%%%%%%%%%%%%%%%%%
% Input files: Please search your DIA data by DIA-NN (*MBR option should be uncheck) 
% and use the "report.tsv" to organize and manupulate the results in a new "csv" file 
% according to the example "data_table" below (imported data in the 1st step)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

close all
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% *Step 1) Importing Data
data_table = readtable('DIANNv18_rawoutput_B57_37C.csv');

rawdata = table2array(data_table(:,[3:8]));
rawrep = table2array(data_table(:,1));
rawseq = table2array(data_table(:,2));

rep_ref = {'37C_r1','37C_r2','37C_r3','42C_r1','42C_r2','42C_r3',...
    '46C_r1','46C_r2','46C_r3','50C_r1','50C_r2','50C_r3','54C_r1','54C_r2','54C_r3',...
    '58C_r1','58C_r2','58C_r3','63C_r1','63C_r2','63C_r3','68C_r1','68C_r2','68C_r3',...
    '73C_r1','73C_r2','73C_r3'}';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% *Step 2) Importing target list of the peptide sequences
seq_tar = csvimport('DIANNv18_seqtarget_B57_37C.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% *Step 3) Selecting the quantification level for profiling as follows:

% Note: 'y' denotes YES and 'n' denotes NO
% Here we used MS1 peak area for quantitative peptidomics and
% thermostability profiling

% 1) Precursor quantity (calculated by DIA-NN) 
prec_quant = 'n';
% 2) MS1 level peak area
MS1_PA = 'y';
% 3) MS2 level peak areas (intensities)
MS2_Frag = 'n';
% 4) Total peak areas (MS1+MS2), i.e., sum of XICs
Total_raw = 'n';


if prec_quant=='y'
    data = rawdata(:,1);
elseif MS1_PA=='y'
    data = rawdata(:,2);
elseif MS2_Frag=='y'
    data = rawdata(:,3);
elseif Total_raw=='y'
    data = rawdata(:,4);
else
    data = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% *Step 4) Manipulating data to preprocess and generate a "data matrix" 
%%% of the quantified precursors across data points:

for dp = 1:1:size(seq_tar,1)
refseq = seq_tar(dp,1);
idx = find(ismember(rawseq,refseq));
reppick = rawrep(idx,1);
datatar = data(idx,1);
quant = zeros(1,size(rep_ref,1));
for rp = 1:1:size(rep_ref,1)
reptar = rep_ref(rp,1);
idr = find(ismember(reppick,reptar));
chimspec = datatar(idr,1);
if size(idr,1)>1
    prec = sum(chimspec(:,1));
elseif size(idr,1)==1 
    prec = chimspec;
else
    prec = 0;
end
quant(1,rp) = prec;
idr = [];
end
data_quant(dp,:) = quant;
idx = []; 
dp
end

seq = seq_tar;
dataini = data_quant;
seg = size(dataini,2)/3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% *Step 5) Isolating iRT peptides' data for the normalization to iRTs:

iRTs_norm='y';

iRT_seq = {'LGGNEQVTR','GAGSSEPVTGLDAK','VEATFGVDESNAK','YILAGVENSK',...
'TPVISGGPYEYR','TPVITGAPYEYR','DGLDAASYYAPVR','ADVTPADFSEWSK',...
'GTFIIDPGGVIR','GTFIIDPAAVIR','LFLQFGAQGSPFLK'}';

if iRTs_norm=='y'
  
       
for di=1:1:size(iRT_seq,1)
ref_irt = iRT_seq{di,1};
idirt = find(ismember(seq, ref_irt));
DiRT(di,:) = dataini(idirt,:);
idirt = [];
end 


%%%%%%% iRTs normalization %%%%%%%

iRT_avr = mean(DiRT,2);
DiRT_norm = DiRT./iRT_avr;
DiRT_avr = mean(DiRT_norm,1);

dataini = dataini./DiRT_avr;

end

data = dataini;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% *Step 6) Rolling the data points to smooth the profile and manipulating 
%%% outliers with constant values

for i = 1:1:seg
minc = (i*3)-2;
maxc = (i*3);
a = data(:,minc:maxc);    
 
for k=1:1:size(a,1)
zc(k) = sum(a(k,:)==0);
end

for j=1:1:size(a,1)
if zc(j)==1
ac = a(j,:);    
idz = find(ac==0);
ac(:,idz)=[];

ma = median(ac);
mn = min(ac);
mx = max(ac);

elseif zc(j)==0
sc = sort(a(j,:),'descend');
ac = sc(:,[1:3]);  

ma = median(ac);
mn = min(ac);
mx = max(ac);

else
ma = 0;
mn = 0;
mx = 0;
end
mia(j) = ma;
mna(j) = mn;
mxa(j) = mx;
end

outmedian(:,i) = mia';
outmin(:,i) = mna';
outmax(:,i) = mxa';
end

for rm = 1:1:size(outmedian,1)
as = outmedian(rm,:);
asmn = outmin(rm,:);
asmx = outmax(rm,:);
for r = 1:1:(size(as,2)-1)
asr = (as(r)+as(r+1))/2;
as_rol(r,1) = asr;
asr = [];

asrmn = (asmn(r)+asmn(r+1))/2;
asmn_rol(r,1) = asrmn;
asrmn = [];

asrmx = (asmx(r)+asmx(r+1))/2;
asmx_rol(r,1) = asrmx;
asrmx = [];
end

if as(size(as,2))>as_rol(end)
lp = 0; %last data point at 73C    
else
lp = as(size(as,2));
end

if asmn(size(asmn,2))>asmn_rol(end)
lpmn = 0; %last data point at 73C    
else
lpmn = asmn(size(asmn,2));
end

if asmx(size(asmx,2))>asmx_rol(end)
lpmx = 0; %last data point at 73C    
else
lpmx = asmx(size(asmx,2));
end

as_rol = [as_rol;lp]';
medrol(rm,:) = as_rol;
asmn_rol = [asmn_rol;lpmn]';
minrol(rm,:) = asmn_rol;
asmx_rol = [asmx_rol;lpmx]';
maxrol(rm,:) = asmx_rol;

as = [];
as_rol = [];
asmn = [];
asmn_rol = [];
asmx = [];
asmx_rol = [];
end


ref = medrol(:,1);
refmn = minrol(:,1);
refmx = maxrol(:,1);
D = medrol(:,1:end);
Dmn = minrol(:,1:end);
Dmx = maxrol(:,1:end);


tol = 1.0;

for ic = 1:1:size(D,1)   
if D(ic,7)>(tol*D(ic,6))
D(ic,[7]) = mean([D(ic,[6]),D(ic,[7])]);
Dmn(ic,[7]) = mean([Dmn(ic,[6]),D(ic,[7])]);
Dmx(ic,[7]) = mean([Dmx(ic,[6]),D(ic,[7])]);
end
end

for ic = 1:1:size(D,1)   
if D(ic,8)>D(ic,7)
D(ic,[8]) = 0;
Dmn(ic,[8]) = 0;
Dmx(ic,[8]) = 0;
end
end

for ic = 1:1:size(D,1)   
if D(ic,9)>D(ic,8)
D(ic,[9]) = 0;
Dmn(ic,[9]) = 0;
Dmx(ic,[9]) = 0;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% *Step 7) Removing precursors with missing values at the any first datapoints
% (37C to 46C)

for h=1:1:size(D,1)
zrid0(h) = sum(medrol(h,1:2)==0);
end

idsel0 = find(zrid0==0);
Dsel0 = D(idsel0,:);
Dmn0 = Dmn(idsel0,:);
Dmx0 = Dmx(idsel0,:);
ref0 = ref(idsel0,:);
data0 = data(idsel0,:);

% Normalizing precursors intensities to 37C data point
Dnrm = Dsel0./(Dsel0(:,1));
Dmn_nrm = Dmn0./(Dsel0(:,1));
Dmx_nrm = Dmx0./(Dsel0(:,1));
seq0 = seq(idsel0,1);
Kb = Dsel0./ref0;
Kd = 1-Kb;

outmedian = outmedian(idsel0,:);


%%% *Step 8) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Sigmoid curve fitting to estimate Tm value (IC50) %%%%
%%%% for thermostability profiling of immunopeptidomes %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Definng a list, showing apllied temperatures in the experiments
tmp = [37,42,46,50,54,58,63,68,73];

Y = Dnrm;
Ymn = Dmn_nrm;
Ymx = Dmx_nrm;
alpha = 1.0;
[m n] = find(Y(:,:)>alpha);
outl = unique(m,'rows');

x1 = tmp';

Yfit = zeros(size(Y,1),100);
PCC = zeros(size(Y,1),1);

Yt = Y;

xx = linspace(min(x1),max(x1),100);

for f=1:1:size(Yt,1)

y1 = Yt(f,:);

[cf,G,fitout,fo,ft]=L4P(x1,y1');

Tmi(f) = cf.C;
corcof(f) = G.rsquare;
yf = cf(xx);
yfit(f,:)=yf';

clc
f
disp('progress%')
prog = (f/size(Yt,1))*100
end

Tm = zeros(1,size(Yt,1));
ttr = xx;

Tm=Tmi;

Tm = Tm';
PCC = corcof';
seqt = seq0;
Yfit = yfit;

%%% *Step 9)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Filtering valid Tm values based on the expected range and
%%%% removing failed fittings 

ftid = find(Tm>37 & Tm<73);

Tm = Tm(ftid,1);
data_sel = data0(ftid,:);
outmedian = outmedian(ftid,:);

Yfit = Yfit(ftid,:);
PCC = PCC(ftid,1);

Y = Y(ftid,:);
Ymn = Ymn(ftid,:);
Ymx = Ymx(ftid,:);
Kb = Kb(ftid,:);
Kd = Kd(ftid,:);

seqt = seqt(ftid,1);

%%% *Step 10)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%***** Evaluation/Validation of the fitting quality *****%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Percentage of the fitted curves with a good correlation with 
% experimental data points, R >= 0.75  

fit_cf = 0.75; % R >= 0.75
   
idfit = find(PCC>=fit_cf);
disp('Percentage% of the fitted curves with a good correlation with experimental data points')
fit_per = (size(idfit,1)/size(PCC,1))*100

%%% *Step 11) Basic visualization to check results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting a histogram to show the distribution for the calculated Tm
% values

figure(1);
histogram(Tm,[37:2:73]); % "Tm" contains the refined list of valid Tm values
ylabel('# profiled peptides');
xlabel('Temperature (C)');
xticks([37:2:73]);

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% *Step 12) Drawing Denaturation Profiles with the correspong Tm values
%%% for individual HLA peptides
%%% e.g., HLA-B*57:01-bound 9-mer peptide; 'ASLQNRIDW'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% *NOTE: If you need to check the denaturation  
%%% profile for another sequence, check the 'seqt' vector from the Workspace
%%% and change the sequence of "intseq" - at line 418 - to visualize another peptide 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grc = [0.93,0.93,0.93];
gray = [0.5,0.5,0.5];

intseq = {'ASLQNRIDW'};

fity = Yfit;
sq = find(ismember(seqt, intseq));

rest = abs(xx-Tm(sq,1));
[mt nt] = min(rest);


ylin = fity(sq,nt);
ymin = Y(sq,:)-Ymn(sq,:);
ymax = Ymx(sq,:)-Y(sq,:);

figure(2);
patch([tmp fliplr(tmp)].',[Ymn(sq,:) fliplr(Ymx(sq,:))].',grc,'EdgeColor','none');
hold on;
errorbar(tmp,Y(sq,:),ymin,ymax,'Marker','o','MarkerEdgeColor','b','LineStyle','none');
hold on;
plot(xx,fity(sq,:),'color','g');
hold on;
plot(Tm(sq,1),ylin,'*','MarkerFaceColor','m','MarkerEdgeColor','m');
hold on;
line([min(xlim()),Tm(sq,1)], [ylin,ylin], 'LineWidth', 0.3, 'Color', 'k', 'LineStyle',':');
hold on;
line([Tm(sq,1),Tm(sq,1)], [min(ylim()),ylin], 'LineWidth', 0.3, 'Color', 'k', 'LineStyle',':');

hold on; 

title({[seqt{sq,1}]
       ['Tm = ' num2str(round(Tm(sq,1),2))]
       });
ylabel('Normalized Peak Area');
xlabel('Temperature (C)');
xticks([37 42 46 50 54 58 63 68 73]);
xticklabels({'37C','42C','46C','50C','54C','58C','63C','68C','73C'});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Export important results, e.g., Tm values, profiled peptide sequences, and fitted curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% User can save "Tm", "seqt", and "Yfit" as the *csv files %%% 
%%% Example (you can change the file name and directory (file path) and
%%% uncomment the following lines to run them and save the main outputs as 'csv'!


% label_param = {'Tm'};
% T = array2table(Tm,'RowNames',seqt,'VariableNames',label_param);
% 
% writetable(T,'D:\Projects_results\C1R_A2B57\Tm_C1R-B5701_FusionData.csv','Delimiter',',');
% writecell(seqt,'D:\Projects_results\C1R_A2B57\Sequences_C1R-B5701_FusionData.csv','Delimiter',',');
% writematrix(Yfit,'D:\Projects_results\C1R_A2B57\Yfit_curves_C1R-B5701_FusionData.csv','Delimiter',',');
