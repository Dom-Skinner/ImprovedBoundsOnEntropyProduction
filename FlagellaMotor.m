clear
load('MTB24_10mM')
file_struct = f;
load('MTB24_85mM')
file_struct = [file_struct,f];
load('MTB32_85mM')
file_struct = [file_struct,f];


%% Calculate sig1 and sig2 from the 
for idx = 1:length(file_struct)
    x_traj = file_struct(idx).x_traj;  
    y_traj = file_struct(idx).y_traj;

    x_mean = mean(x_traj);
    y_mean = mean(y_traj);
    theta = atan2(y_traj-y_mean,x_traj - x_mean);
    theta = theta_correct(theta);
    increment = 1e-04; % same for all experiments.
    t = increment*length(theta);
    
    [sig1,sig2] = entropy_sweep(theta,1:15,t);

    file_struct(idx).sig1 = max(sig1);
    file_struct(idx).sig2 = max(sig2);
end


%% Make some plots
keep_MTB24_10mM = extract_type(file_struct,"MTB24_10mM");
sig1_MTB24_10mM = [file_struct(keep_MTB24_10mM).sig1];
sig2_MTB24_10mM = [file_struct(keep_MTB24_10mM).sig2];

keep_MTB24_85mM = extract_type(file_struct,"MTB24_85mM");
sig1_MTB24_85mM = [file_struct(keep_MTB24_85mM).sig1];
sig2_MTB24_85mM = [file_struct(keep_MTB24_85mM).sig2];

keep_MTB32_85mM = extract_type(file_struct,"MTB32_85mM");
sig1_MTB32_85mM = [file_struct(keep_MTB32_85mM).sig1];
sig2_MTB32_85mM = [file_struct(keep_MTB32_85mM).sig2];

%%{
x1 = [sig1_MTB24_10mM';sig1_MTB24_85mM'; sig1_MTB32_85mM'];
x2 = [sig2_MTB24_10mM';sig2_MTB24_85mM'; sig2_MTB32_85mM'];


g1 = repmat({'MTB24 10mM'},length(sig1_MTB24_10mM),1);
g2 = repmat({'MTB24 85mM'},length(sig1_MTB24_85mM),1);
g3 = repmat({'MTB32 85mM'},length(sig1_MTB32_85mM),1);
g4 = repmat({'MTB24 10mM s2'},length(sig1_MTB24_10mM),1);
g5 = repmat({'MTB24 85mM s2'},length(sig1_MTB24_85mM),1);
g6 = repmat({'MTB32 85mM s2'},length(sig1_MTB32_85mM),1);

g = [g1; g2;g3;g4;g5;g6];
boxplot([x1,x2],g)
ylabel('Rate of entropy production')

function theta = theta_correct(theta)
    for i = 2:length(theta)
        if abs(theta(i) - theta(i-1)) > 1.5*pi
            theta(i) = theta(i) - 2*pi*round( (theta(i)-theta(i-1))/2/pi);
        end
    end
end

function keep = extract_type(file_struct,type)
    keep = true(length(file_struct),1);
    for i = 1:length(file_struct)
        keep(i) = (file_struct(i).type == type);
    end
end