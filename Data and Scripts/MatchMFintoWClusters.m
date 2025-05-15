%%%
%%% This is to match all the individual muscle field into synergy clusters
%%% Input:
%%%     "cluster_ms" contains the clusters of muscle synergies
%%%     "MPI_total_new", the sub-structure MPI_total_new(X).SS_T.mpi
%%%     contains the normalized muscle field of units.
%%% 
%%%
%%% This script concerns only the vectors from single stimulation at 
%%% threshold level (SST) condition. It could be specified in the 
%%% "conditions".
%%% 
%%%
%%% Created by Borong @CUHK, 8th Oct 2022

clc
clear
close all

%  ------------------------------------------------------------------------
%  Stack muscle fields
%  ------------------------------------------------------------------------
emg_names_new = [{'TA'},{'GC'},{'BF'},{'VL'},{'VM'},{'GM'}];
MPI_data = matfile('MPI_total_new.mat');
MPI_total_new = MPI_data.MPI_total_new;
cluster_ms_data = matfile('cluster_ms.mat');
cluster_ms = cluster_ms_data.cluster_ms;
preference_data = matfile('preference_pair.mat');
preference_pair = preference_data.preference_pair;
clear MPI_data cluster_ms_data preference_data
condition = [{'SS_T'}];
MF_total = [];
MF_index = [];
for i = 1:length(MPI_total_new)
    for c = 1:length(condition)
        if ~isempty(MPI_total_new(i).(condition{c}))
            data = MPI_total_new(i).(condition{c}).mpi;
            for d = 1:size(data,2)
                mf_temp(:,d) = data(:,d);
                % Normalize MF vectors inside "data" structure, 'mf_temp'
                % is column vectors
                mf_temp(:,d) = mf_temp(:,d) ./ sqrt(mf_temp(:,d)'* mf_temp(:,d));
            end
            MF_total = [MF_total ; transpose(mf_temp)];
            MF_index = [MF_index ; i*ones(size(mf_temp,2),1)];
            clear mf_temp d data
        end
    end
    clear c
end
clear i


%  ------------------------------------------------------------------------
%  Match muscle field into clusters of muscle synergies 
%  Input: "MF_total" with dimensions of number of muscle fields and number 
%         of muscles.
%         "cluster_ms" contains the clusters of muscle synergies.
%  Output: "id_MFtoWClus" is a column vector where each row corresponds to
%         the matched W clusters.
%  ------------------------------------------------------------------------

for m = 1:size(MF_total,1)
    data_mf = MF_total(m,:);
    if ~isnan(data_mf)
        mf_ind = MF_index(m);
        % Match the mth muscle field into its mth preferred synergy cluster
        for c = 1:length(cluster_ms)
            subject = cluster_ms{c}.subject;
            data_ms = cluster_ms{c}.data;
            if ~isempty(find(mf_ind==subject, 1))
                id(m) = find(mf_ind==subject, 1);
                sp_temp(c) = dot(data_mf,data_ms(id(m),:));
            end
            clear data_ms subject
        end
        [~,id_clus(m)] = max(sp_temp);
        clear c mf_ind sp_temp


%         % Calculate the similarity of one muscle field with synergies of 
%         % each cluster of synergies.
%         for c = 1:length(cluster_ms) 
%             data_ms = cluster_ms{c}.data;
%             for d = 1:size(data_ms,1) % "data_ms" is a row column
%                 % Compare the similarity with each of the synergy within
%                 % the cth cluster of synergies.
%                 sp_temp(d) = dot(data_mf,data_ms(d,:));
%             end
%             sp_mf_WClus(c) = mean(sp_temp);
%             clear sp_temp d
%         end
%         [~,id_clus(m)] = max(sp_mf_WClus);
%         clear sp_mf_WClus data_ms c
    end
end
clear data_mf m


%% ------------------------------------------------------------------------
%  Plot cluster of synergies with their matched muscle fields.
%  ------------------------------------------------------------------------

for i = 1:length(cluster_ms)
    W = cluster_ms{i}.data'; % "W" is a column vector
    MF = MF_total(id_clus==i,:)'; % "MF" is a column vector
    sub_mf_temp = MF_index(id_clus==i);
    NMus = size(W,1);
    H(i) = figure,
    for j=1:size(W,2)+size(MF,2)
        fighandle{j} = subplot(1,size(W,2)+size(MF,2),j,'align');    
        set(gca,'ytick',1:NMus,'ylim',[0 NMus+1],'xlim',[0 1],...
            'xtick',[0 1],'xticklabel',[],'ygrid','on','box','off')
        if le(j,size(W,2))
            sub_ms = cluster_ms{i}.subject(j);
            barh(flipud(W(:,j)),0.4,'FaceColor',[0 0 0.5]);
            title(num2str(sub_ms));
        else
            sub_mf = sub_mf_temp(j-size(W,2));
            barh(flipud(MF(:,j-size(W,2))),0.4,'FaceColor',[0.5 0 0]);
            title(num2str(sub_mf));
        end    
        if isequal(j,1)
            set(gca,'yticklabel',fliplr(emg_names_new));
            xticklabels([]);
        else
            set(gca,'yticklabel',[]);
            xticklabels([]);
        end
        set(gca,'FontSize',20);
    end
end











