function [ sta ] = compute_sta( stim, rho, para )
        

        %COMPUTE_STA Calculates the spike-triggered average for a neuron that
        %            is driven by a stimulus defined in stim. The spike-
        %            triggered average is computed over num_timesteps timesteps.
            
        % Should add an option of whether having 'para.pre' and 'para.post'
        % input
        
        if para.pre_post
            pre_timesteps = para.pre_timesteps;
            post_timesteps = para.post_timesteps;
            
            num_timesteps = para.num_timesteps;
            sta = zeros(num_timesteps, 1);
            
        %    spike_times = find(rho(num_timesteps+1:end)) + num_timesteps;
        %    spike_times = find(rho(1:end-num_timesteps));
            
            spike_times = find(rho(pre_timesteps:end-post_timesteps));

            % Compute the spike-triggered average of the spikes found using the
            % find command. To do this, compute the average of all of the vectors
            % starting from a spike and ending at the time after 15ms of
            % the event (inclusive). 
            % Each of these vectors defines a list of
            % samples that is contained within a window of 15 ms after the each
            % spike. 
            

            [~,I] = find(rho==1);
            stim_new = zeros(length(I),num_timesteps); % all the indices of SPT
            for i = 1:length(I)
                if ge(I(i),pre_timesteps) && le(I(i),size(stim,2)-post_timesteps)
                    stim_new(i,:) = stim(I(i)-pre_timesteps:I(i)+post_timesteps-1)';
                end
            end

            sta = mean(stim_new,1);
            
        else
            
            num_timesteps = para.num_timesteps;
            sta = zeros(num_timesteps, 1);
            spike_times = find(rho(1:end-num_timesteps));
            
            
            [~,I] = find(rho==1);
            stim_new = zeros(length(I),num_timesteps); % all the indices of SpTA
            for i = 1:length(I)
                if le(I(i),size(stim,2)-num_timesteps)
                    stim_new(i,:) = stim(I(i):I(i)+num_timesteps-1)';
                end
            end

            sta = mean(stim_new,1);
        end



            



end
