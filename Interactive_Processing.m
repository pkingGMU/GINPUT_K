clear all
close all
clc

%directory = '/Users/kpore/Documents/PhD Kinesiology/Research/Raffegeau/SwitzerFellowship/Data';
% directory = '/Users/kpore/Documents/PhD Kinesiology/Research/Raffegeau/PDTVR/HMDData/1minPaper';
% directory = '/Users/kpore/Documents/PhD Kinesiology/Research/Raffegeau/Kelly_PilotDiss/NewHeadset';
% directory = '/Users/sarah/Box/data analyst contract work/raffegeau gmu/DT VR in YA/data/HMD';
directory = strcat(pwd, '/Data');

fileList = dir(directory);
fileList([fileList.isdir],:)= [];
fileList = fileList(contains({fileList.name}, '.csv'));
% fileList=fileList(5:end, :);

% for gg = 1:24;
for gg = 1:2
% for gg = 51
clear RHS LHS LTO RTO rmintab rmaxtab lmintab lmaxtab ltoeloc ltoemin rtoeloc rtoemin allevents

    filename = fileList(gg).name;
    % 2 trials have no step data (just standing still)
%     if ~strcmp(filename,'_SUB_DT VR 009_2_HIGH_ST_(12-02-21 11-43-04).csv') & ~strcmp(filename,'_SUB_DT VR 018_2_High_DT_(18-06-21 15-43-18).csv') & ~strcmp(filename,'_SUB_DTVR013_5_HIGH_DT_(08-04-21 16-04-52).csv') & ~strcmp(filename,'_SUB_____(10-02-21 18-25-17).csv')
        %% Set up the Import Options and import the data
        opts = delimitedTextImportOptions("NumVariables", 16);
        opts.ExtraColumnsRule='ignore';

        % Specify range and delimiter
        opts.DataLines = [1, 2];
        opts.Delimiter = ",";

        %Import two data tables
        headerInfo = readtable(strcat(directory,"/",filename), opts);
        clear opts;

        headerInfo.Properties.VariableNames = ["Subject","Trial", "Task", "Elevation", "Corner_BottomLeft_PositionX", "Corner_BottomLeft_PositioY", "Corner_BottomLeft_PositionZ",...
            "Corner_BottomRight_PositionX", "Corner_BottomRight_PositioY", "Corner_BottomRight_PositionZ", "Corner_TopRight_PositionX","Corner_TopRight_PositioY", "Corner_TopRight_PositionZ", "Corner_TopLeft_PositionX","Corner_TopLeft_PositioY","Corner_TopLeft_PositionZ"];
        %     if gg < 450
        %
        %     else
        %         headerInfo.Properties.VariableNames = ["Subject","Trial", "Task", "Elevation", "Corner_BottomLeft_PositionX", "Corner_BottomLeft_PositioY", "Corner_BottomLeft_PositionZ",...
        %             "Corner_BottomRight_PositionX", "Corner_BottomRight_PositioY", "Corner_BottomRight_PositionZ", "Corner_TopRight_PositionX","Corner_TopRight_PositioY", "Corner_TopRight_PositionZ", "Corner_TopLeft_PositionX","Corner_TopLeft_PositioY","Corner_TopLeft_PositionZ",...
        %             "My_thoughts_are", "My_body_feels",	"I_am_feeling",	"Mental_effort"];
        %     end

        headerInfo = headerInfo(2,:);
        data = readtable(strcat(directory,"/",filename), 'HeaderLines', 2);
        %% define plank coordinates
        %plank = [x,y,z]
        BottomLeft = str2double([headerInfo.Corner_BottomLeft_PositionX, headerInfo.Corner_BottomLeft_PositioY,...
            headerInfo.Corner_BottomLeft_PositionZ]);

        BottomRight= str2double([headerInfo.Corner_BottomRight_PositionX, headerInfo.Corner_BottomRight_PositioY,...
            headerInfo.Corner_BottomRight_PositionZ]);

        TopLeft = str2double([headerInfo.Corner_TopLeft_PositionX, headerInfo.Corner_TopLeft_PositioY,...
            headerInfo.Corner_TopLeft_PositionZ]);

        TopRight= str2double([headerInfo.Corner_TopRight_PositionX, headerInfo.Corner_TopRight_PositioY,...
            headerInfo.Corner_TopRight_PositionZ]);

        %%%Getting trajectories for sensors in linear matrix

        %Left foot sensor
        LFootlin = data(:,21:23);
%         LFootlin = data(:,27:29);
        %LFootrot = [data.Foot_1RotationX, data.Foot_1RotationY, data.Foot_1RotationZ];
        %Right foot sensor
        RFootlin = data(:,27:29);
%         RFootlin = data(:,21:23);
        %RFootrot = [data.Foot_2RotationX, data.Foot_2RotationY, data.Foot_2RotationZ];
        %HMD sensor
        HMDlin = data(:,3:5);

        %HMDrot = [data.HeadSetRotationX, data.HeadSetRotationY,...
        %   data.HeadSetRotationZ];
        LFootlin = table2array(LFootlin);
        RFootlin = table2array(RFootlin);
        HMDlin = table2array(HMDlin);

        %create average vectors from existing coordinates
        origin = mean([BottomLeft; BottomRight],1);
        endpoint = mean([TopLeft; TopRight],1);

        %set direction of vectors
        x = (endpoint - origin);
        z = (BottomRight - BottomLeft);
        y = cross(z,x);

        %normalize z using x because its most reliable, using x,y to orient vector in
        %proper direction
        z = cross (x,y);


        %logic statement to check the vectors are orthogonal
        dot(z,y) == 0 & dot(z,x) == 0;

        %determine magnitude of vectors
        i = x ./ norm(x);
        j = y ./ norm(y);
        k = z ./ norm(z);

        G2A = [i; j; k;];
        %A2G = G2A';

        xDir = [origin; i+origin];
        yDir = [origin; j+origin];
        zDir = [origin; k+origin];

        %%Redefine plank positions
        %subtract origin value from all marker data
        %transposed first for matrix math and transposed again for graphing
        BLeft = [G2A*[BottomLeft-origin]']';
        BRight = [G2A*[BottomRight-origin]']';
        TLeft = [G2A*[TopLeft-origin]']';
        TRight = [G2A*[TopRight-origin]']';

        %plot coordinate system back in local frame, not global
        i2 = [G2A*i']';
        j2 = [G2A*j']';
        k2 = [G2A*k']';

        %bringing back to local 000 origin not the newly defined origin for plot
        xDir = [[0,0,0]; i2];
        yDir = [[0,0,0]; j2];
        zDir = [[0,0,0]; k2];

        %% Creating corrected trajectories
        %linear
        lfoot = [G2A*[LFootlin-origin]']';
        rfoot = [G2A*[RFootlin-origin]']';
        HMD = [G2A*[HMDlin-origin]']';
        %%flipping the axes for so that vertical is positive
        HMD = [HMD.*[1, -1, 1]];
        lfoot = [lfoot.*[1, -1, 1]];
        rfoot = [rfoot.*[1, -1, 1]];

        %% interpolate signal for equal time sample
        %%HMD
        HMD_time = timeseries(HMD, data.ElapsedTime);
        new_interp_time = [0:0.01:data.ElapsedTime(end,:)];
        HMD_int = resample(HMD_time, new_interp_time);

        %get rid of NAN in data set (first few rows) replace with 0
        HMD_int.Data(isnan(HMD_int.Data)) = 0;

        %changes time series to a matrix for later use
        HMD_corrected = HMD_int.Data;

        %%rfoot
        rfoot_time = timeseries(rfoot, data.ElapsedTime);
        rfoot_int = resample(rfoot_time, new_interp_time);

        %get rid of NAN in data set (first few rows) replace with 0
        rfoot_int.Data(isnan(rfoot_int.Data)) = 0;

        %changes time series to a matrix for later use
        rfoot_corrected = rfoot_int.Data;
        rfoot_final = rfoot_corrected;
        %for gap filling
        rfoot_smooth8 = movmean(rfoot_corrected, 8);

        %%lfoot
        lfoot_time = timeseries(lfoot, data.ElapsedTime);
        lfoot_int = resample(lfoot_time, new_interp_time);

        %get rid of NAN in data set (first few rows) replace with 0
        lfoot_int.Data(isnan(lfoot_int.Data)) = 0;

        %changes time series to a matrix for later use
        lfoot_corrected = lfoot_int.Data;
        lfoot_final = lfoot_corrected;

        %for gap filling
        lfoot_smooth8 = movmean(lfoot_corrected, 8);

        %     %%plots for a job talk
        %     figure();
        %     plot(BLeft, BRight, TLeft, TRight)
        %     hold on
        %     plot(origin)
        %     %draw(rfoot_corrected)
        %     %draw(lfoot_corrected)
        %     %draw(HMD_corrected)
        %     hold off



        %% Gap Filling
        %%this is not always needed, especially not in short trials
        %%find parts of signal that are ruined by poor data tracking
        fs = 100;
        %%creating an error threshold for the values 9/15/21
        %     errorThresholdL = -origin(:,2);%mean(lfoot_smooth8(fs*0.5:fs*3.5,2))*0.8% - *std(lfoot_corrected(fs*0.5:fs*3.5,2))
        %     errorThresholdR = -origin(:,2);%mean(rfoot_smooth8(fs*0.5:fs*3.5,2))*0.8% - 8*std(lfoot_corrected(fs*0.5:fs*3.5,2))
        %%if high add 15m to error threshold..
        if any(lfoot_smooth8(:,2)>15)
            errorThresholdL= -origin(:,2)+15;
        else
            errorThresholdL= -origin(:,2);
        end

        if any(rfoot_smooth8(:,2)>15)
            errorThresholdR= -origin(:,2)+15;
        else
            errorThresholdR= -origin(:,2);
        end

        %%identifying time points where the data is beyond said error threshold
        clear s s2
        %%left foot
        k = 1;
        vrTime = lfoot_int.Time;
        [badL] = find(lfoot_smooth8(:,2)<errorThresholdL);
        if isempty(badL)
            lgapFill = false;
        else
            gap = find(diff(badL)>k*1);

            if isempty(gap)
                gap = length(badL);
            end

            badL_wGap = [];
            for n=1:length(gap)
                if n==1
                    badL_wGap(n,:) = [max([badL(1)-k; 1]) badL(gap(n))+k];
                else
                    badL_wGap(n,:) = [badL(gap(n-1)+1)-k badL(gap(n))+k];
                end
            end

            badL_complete = [];
            for n=1:length(badL_wGap(:,1))
                badL_complete = [badL_complete;[badL_wGap(n,1):badL_wGap(n,2)]'];
            end

            lfoot_interpolated = [vrTime, lfoot_smooth8];
            lfoot_interpolated(badL_complete,:) = [];


            for n=1:3
                s(:,n) = interp1(lfoot_interpolated(:,1), lfoot_interpolated(:,n+1), vrTime(badL_complete), 'l');
            end
            lfoot_final(badL_complete,:) = s;
            lgapFill = true;
        end

        %%rfoot
        vrTime = rfoot_int.Time;
        [badR] = find(rfoot_smooth8(:,2)<errorThresholdR);
        if isempty(badR)
            rgapFill = false;
        else
            gap = find(diff(badR)>k*1);

            if isempty(gap)
                gap = length(badR);
            end

            badR_wGap = [];
            for n=1:length(gap)
                if n==1
                    badR_wGap(n,:) = [max([badR(1)-k; 1]) badR(gap(n))+k];
                else
                    badR_wGap(n,:) = [badR(gap(n-1)+1)-k badR(gap(n))+k];
                end
            end

            badR_complete = [];
            for n=1:length(badR_wGap(:,1))
                badR_complete = [badR_complete;[badR_wGap(n,1):badR_wGap(n,2)]'];
            end

            rfoot_interpolated = [vrTime, rfoot_smooth8];
            rfoot_interpolated(badR_complete,:) = [];

            for n=1:3
                s2(:,n) = interp1(rfoot_interpolated(:,1), rfoot_interpolated(:,n+1), vrTime(badR_complete), 'l');
            end
            %to compare gap-filled trials with same trial, remove
            rfoot_final(badR_complete,:) = s2;
            rgapFill = true;
        end
        clear badR_complete badL_complete


        %     rfoot_final(~any(lfoot_final,2),:) = [];
        %     lfoot_final(~any(rfoot_final,2),:) = [];
        %     HMD_corrected(~any(rfoot_final,2),:) = [];
        %



        %%ChatGPT solution%%
        %woot!
        % Find the rows in rfoot_final that have at least one non-zero element
        non_empty_rows_rfoot_final = any(rfoot_final,2);

        % Find the rows in lfoot_final that have at least one non-zero element
        non_empty_rows_lfoot_final = any(lfoot_final,2);

        % Find the rows in rfoot_final that have at least one non-zero element
        non_empty_rows_HMD_corrected = any(HMD_corrected,2);

        % Find the common rows that have at least one non-zero element in all three matrices
        common_non_empty_rows = non_empty_rows_rfoot_final & non_empty_rows_lfoot_final & non_empty_rows_HMD_corrected;

        % Extract the common rows from all three matrices
        rfoot_final = rfoot_final(common_non_empty_rows,:);
        lfoot_final = lfoot_final(common_non_empty_rows,:);
        HMD_corrected = HMD_corrected(common_non_empty_rows,:);

        %filtering parameters
        fc = 6;
        [b,a] = butter(2,fc/(fs/2));

        %filter corrected trajectory
        HMD_filt= filtfilt(b,a,HMD_corrected);

        %%rfoot
        rfoot_filt= filtfilt(b,a,rfoot_final);
        %%lfoot
        lfoot_filt= filtfilt(b,a,lfoot_final);

        %%use velocity of HMD loop calculator here
        velHead = zeros(length(HMD_filt),3);
        velRFoot = zeros(length(rfoot_filt),3);
        velLFoot = zeros(length(lfoot_filt),3);

        for n=2:length(velHead)-1
            velHead(n,:) = (HMD_filt(n+1,:) - HMD_filt(n-1,:)) ./(new_interp_time(n+1) - new_interp_time(n-1));
        end

        for n=2:length(velRFoot)-1
            velRFoot(n,:) = (rfoot_filt(n+1,:) - rfoot_filt(n-1,:)) ./(new_interp_time(n+1) - new_interp_time(n-1));
        end

        for n=2:length(velLFoot)-1
            velLFoot(n,:) = (lfoot_filt(n+1,:) - lfoot_filt(n-1,:)) ./(new_interp_time(n+1) - new_interp_time(n-1));
        end

        %%smooth the filtered data to calculate still period better
        lfoot_velocitysmooth = movmean(vecnorm(velLFoot(2:end,1:2)'), 20);
        rfoot_velocitysmooth = movmean(vecnorm(velRFoot(2:end,1:2)'), 20);
       
        %% creating a distance vector between the HMD and the foot tracker
        rGE_vector = (rfoot_filt(:,1)-HMD_filt(:,1)).*sign(velHead(:,1));
        lGE_vector = (lfoot_filt(:,1)-HMD_filt(:,1)).*sign(velHead(:,1));

        %find peaks that represent heel contact-max distance between HMD and
        %tracker

        [rmaxtab, rmintab] = findpeaks(rGE_vector,'MinPeakProminence',0.1);
        [lmaxtab,lmintab] = findpeaks(lGE_vector,'MinPeakProminence',0.1);

        [rtoemin, rtoeloc] = findpeaks(-rGE_vector, 'MinPeakProminence', 0.1);
        [ltoemin, ltoeloc] = findpeaks(-lGE_vector, 'MinPeakProminence', 0.1);


%% Right Plot Selection 

% Step 3: Plot the original signal and detected peaks
figure;
plot(rGE_vector, 'b');  % Plot the right signal in blue
hold on;
plot(rmintab, rmaxtab, 'mo', 'MarkerFaceColor', 'm');  % Mark right peaks with red circles

plot(lGE_vector, 'Color', [0, 0.5, 0]);  % Plot the left signal in darker green
plot(lmintab, lmaxtab, 'co', 'MarkerFaceColor', 'c');  % Mark left peaks with magenta circles

plot(rGE_vector, 'b');  % Plot the right signal in blue
hold on;
plot(rtoeloc, -rtoemin, 'bo', 'MarkerFaceColor', 'b');  % Mark right peaks with red circles

plot(lGE_vector, 'Color', [0, 0.5, 0]);  % Plot the left signal in darker green
plot(ltoeloc, -ltoemin, 'go', 'MarkerFaceColor', 'g');  % Mark left peaks with magenta circles

title('Detected Peaks (Click on peaks to remove from PINK DOTS)');
xlabel('Sample Index');
ylabel('Signal Amplitude');

% Step 4: Allow user to manually select peaks (click on graph)
disp('Click on peaks you want to remove. Press Enter when done.');
[x_manualr, ~] = ginput;  % User clicks on graph to select peaks


% Step 5: Match manually selected points to the closest peaks
tolerance = 60;  % Increased tolerance to allow for better matching
selected_peaksr = false(size(rmintab));  % Logical array to track selected peaks

for i = 1:length(x_manualr)
    % Find the index of the closest peak (in terms of x-values/indices) to the selected point
    [~, idxr] = min(abs(rmintab - round(x_manualr(i))));
    
    % Check if the selected point is within the tolerance range of a detected peak index
    if abs(rmintab(idxr) - round(x_manualr(i))) <= tolerance
        selected_peaksr(idxr) = true;  % Mark peak as selected
    end
end

% Step 6: Filter the peaks to delete the selected ones
rmintab(selected_peaksr) = [];
rmaxtab(selected_peaksr) = [];


%% Left Plot Selection 

% Step 3: Plot the original signal and detected peaks
figure;
plot(rGE_vector, 'b');  % Plot the right signal in blue
hold on;
plot(rmintab, rmaxtab, 'mo', 'MarkerFaceColor', 'm');  % Mark right peaks with red circles

plot(lGE_vector, 'Color', [0, 0.5, 0]);  % Plot the left signal in darker green
plot(lmintab, lmaxtab, 'co', 'MarkerFaceColor', 'c');  % Mark left peaks with magenta circles

plot(rGE_vector, 'b');  % Plot the right signal in blue
hold on;
plot(rtoeloc, -rtoemin, 'bo', 'MarkerFaceColor', 'b');  % Mark right peaks with red circles

plot(lGE_vector, 'Color', [0, 0.5, 0]);  % Plot the left signal in darker green
plot(ltoeloc, -ltoemin, 'go', 'MarkerFaceColor', 'g');  % Mark left peaks with magenta circles

title('Detected Peaks (Click on peaks to remove from CYAN DOTS)');
xlabel('Sample Index');
ylabel('Signal Amplitude');


% Step 4: Allow user to manually select peaks (click on graph)
disp('Click on peaks you want to remove. Press Enter when done.');
[x_manuall, ~] = ginput;  % User clicks on graph to select peaks


% Step 5: Match manually selected points to the closest peaks
tolerance = 60;  % Increased tolerance to allow for better matching
selected_peaksl = false(size(lmintab));  % Logical array to track selected peaks

for j = 1:length(x_manuall)
    % Find the index of the closest peak (in terms of x-values/indices) to the selected point
    [~, idxl] = min(abs(lmintab - round(x_manuall(j))));
    
    % Check if the selected point is within the tolerance range of a detected peak index
    if abs(lmintab(idxl) - round(x_manuall(j))) <= tolerance
        selected_peaksl(idxl) = true;  % Mark peak as selected
    end
end

% Step 6: Filter the peaks to delete the selected ones
lmintab(selected_peaksl) = [];
lmaxtab(selected_peaksl) = [];


% Step 7: Update the plot with the selected peaks
figure;
plot(lGE_vector, 'b');  % Plot the signal
hold on;
plot(lmintab, lmaxtab, 'ro', 'MarkerFaceColor', 'r');  % Mark selected peaks with green circles
title('Manually Selected Peaks');
hold on
plot(rGE_vector, 'b');  % Plot the signal
hold on;
plot(rmintab, rmaxtab, 'bo', 'MarkerFaceColor', 'b');  % Mark selected peaks with green circles
title('Manually Selected Peaks');
plot(rGE_vector, 'b');  % Plot the right signal in blue
hold on;
plot(rtoeloc, -rtoemin, 'bo', 'MarkerFaceColor', 'b');  % Mark right peaks with red circles

plot(lGE_vector, 'Color', [0, 0.5, 0]);  % Plot the left signal in darker green
plot(ltoeloc, -ltoemin, 'go', 'MarkerFaceColor', 'g');  % Mark left peaks with magenta circles


disp('Selected peaks are now displayed in red.');
disp('Selected peaks are now displayed in blue.');



%% Right Toe Plot Selection 

% Step 3: Plot the original signal and detected peaks
figure;
plot(rGE_vector, 'b');  % Plot the right signal in blue
hold on;
plot(rmintab, rmaxtab, 'mo', 'MarkerFaceColor', 'm');  % Mark right peaks with red circles

plot(lGE_vector, 'Color', [0, 0.5, 0]);  % Plot the left signal in darker green
plot(lmintab, lmaxtab, 'co', 'MarkerFaceColor', 'c');  % Mark left peaks with magenta circles

plot(rGE_vector, 'b');  % Plot the right signal in blue
hold on;
plot(rtoeloc, -rtoemin, 'bo', 'MarkerFaceColor', 'b');  % Mark right peaks with red circles

plot(lGE_vector, 'Color', [0, 0.5, 0]);  % Plot the left signal in darker green
plot(ltoeloc, -ltoemin, 'go', 'MarkerFaceColor', 'g');  % Mark left peaks with magenta circles

title('Detected Peaks (Click on peaks to remove from BLUE DOTS)');
xlabel('Sample Index');
ylabel('Signal Amplitude');

% Step 4: Allow user to manually select peaks (click on graph)
disp('Click on peaks you want to remove. Press Enter when done.');
[x_mantoer, ~] = ginput;  % User clicks on graph to select peaks


% Step 5: Match manually selected points to the closest peaks
tolerance = 60;  % Increased tolerance to allow for better matching
selected_peakstoer = false(size(rtoeloc));  % Logical array to track selected peaks

for i = 1:length(x_mantoer)
    % Find the index of the closest peak (in terms of x-values/indices) to the selected point
    [~, idxr] = min(abs(rtoeloc - round(x_mantoer(i))));
    
    % Check if the selected point is within the tolerance range of a detected peak index
    if abs(rtoeloc(idxr) - round(x_mantoer(i))) <= tolerance
        selected_peakstoer(idxr) = true;  % Mark peak as selected
    end
end

% Step 6: Filter the peaks to delete the selected ones
rtoeloc(selected_peakstoer) = [];
rtoemin(selected_peakstoer) = [];


%% Left Plot Selection 

figure;
plot(rGE_vector, 'b');  % Plot the right signal in blue
hold on;
plot(rmintab, rmaxtab, 'mo', 'MarkerFaceColor', 'm');  % Mark right peaks with red circles

plot(lGE_vector, 'Color', [0, 0.5, 0]);  % Plot the left signal in darker green
plot(lmintab, lmaxtab, 'co', 'MarkerFaceColor', 'c');  % Mark left peaks with magenta circles

plot(rGE_vector, 'b');  % Plot the right signal in blue
hold on;
plot(rtoeloc, -rtoemin, 'bo', 'MarkerFaceColor', 'b');  % Mark right peaks with red circles

plot(lGE_vector, 'Color', [0, 0.5, 0]);  % Plot the left signal in darker green
plot(ltoeloc, -ltoemin, 'go', 'MarkerFaceColor', 'g');  % Mark left peaks with magenta circles


title('Detected Peaks (Click on peaks to remove from GREEN DOTS)');
xlabel('Sample Index');
ylabel('Signal Amplitude');

% Step 4: Allow user to manually select peaks (click on graph)
disp('Click on peaks you want to remove. Press Enter when done.');
[x_mantoel, ~] = ginput;  % User clicks on graph to select peaks


% Step 5: Match manually selected points to the closest peaks
tolerance = 60;  % Increased tolerance to allow for better matching
selected_peakstoel = false(size(ltoeloc));  % Logical array to track selected peaks

for j = 1:length(x_mantoel)
    % Find the index of the closest peak (in terms of x-values/indices) to the selected point
    [~, idxl] = min(abs(ltoeloc - round(x_mantoel(j))));
    
    % Check if the selected point is within the tolerance range of a detected peak index
    if abs(ltoeloc(idxl) - round(x_mantoel(j))) <= tolerance
        selected_peakstoel(idxl) = true;  % Mark peak as selected
    end
end

% Step 6: Filter the peaks to delete the selected ones
ltoeloc(selected_peakstoel) = [];
ltoemin(selected_peakstoel) = [];


% Step 7: Update the plot with the selected peaks
figure;
plot(lGE_vector, 'r');  % Plot the signal
hold on;
plot(ltoeloc, -ltoemin, 'ro', 'MarkerFaceColor', 'r');  % Mark selected peaks with green circles
title('Manually Selected Peaks');
hold on
plot(rGE_vector, 'b');  % Plot the signal
hold on;
plot(rtoeloc, -rtoemin, 'bo', 'MarkerFaceColor', 'b');  % Mark selected peaks with green circles
title('Manually Selected Peaks');
plot(lGE_vector, 'r');  % Plot the signal
hold on;
plot(lmintab, lmaxtab, 'ro', 'MarkerFaceColor', 'r');  % Mark selected peaks with green circles
title('Manually Selected Peaks');
hold on
plot(rGE_vector, 'b');  % Plot the signal
hold on;
plot(rmintab, rmaxtab, 'bo', 'MarkerFaceColor', 'b');  % Mark selected peaks with green circles
title('Manually Selected Peaks');
disp('Selected peaks are now displayed in red.');
disp('Selected peaks are now displayed in blue.');

clearvars -except rmintab rmaxtab ltoeloc rtoeloc ltoemin rtoemin lmintab lmaxtab lfoot_final rfoot_final headerInfo ...
    fileList directory
cd("Saved_Process")
save(strcat(string(headerInfo.Subject), string(headerInfo.Trial)));
cd ..



end