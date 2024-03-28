%%% from pipline source construction in MEGtools
%%% 
exp_dir = '/bcbl/home/home_n-z/sgeng/toolbox/toolboxes/';

% tbpath = '~/Dropbox/toolboxes/'
toolboxes = {'cartographie_motrice' 'fieldtrip-20170911' 'freesurfer'};
for tool = toolboxes   
   addpath(genpath(fullfile(exp_dir,tool{1})))  
end


subject = {  
         'S1';
         'S2';
         'S3'        
         'S4';
         'S5';
         'S6';
        'S7';
        'S8'        
        'S9';
          'S11';
          'S12';
       'S13';
        };
    
   
    
 
    raw_exp = 'G:\Premeg\shuang\Freesurfer';
    meg = 'G:\Premeg\shuang\MEG';
    cd(fullfile(raw_exp));
    subjects_dir = fullfile(raw_exp);

    raw_dir =  fullfile(raw_exp,'raw');
    exp = 'trans';




    %% ADD PATHS
    % add your fieldtrip

    ft_defaults 
    % since demean must be yes. 

%% check automaticly which step is not be done.
    for n_sub = 1:length(subject)
        cfg=[]
        cfg.subjects_dir = subjects_dir;
        cfg.subject = subject{n_sub};
        cfg.exp_name = exp;
        cfg.meg_file = fullfile(raw_dir,[subject{n_sub} '_raw_tsss.fif']);
        cfg.dicom = '';  % what is this 
        cfg.space = 5;   % is this the grid size?
        CM_mne_fwd_sol(cfg);
         CM_mne_fwd_sol_from_MNI(cfg);
    end

    


%% Load head info if it is neccessary, normally dont need 
%     load('/bcbl/home/home_n-z/sgeng/projects/naming/Naming_controls/filedtrip_data/sub1/meg/01_Euskera_Objects_tsss_mc.mat')
%     grads = find(strcmp(data.grad.chantype,'megplanar'));
%     magns = find(strcmp(data.grad.chantype,'megmag'));


    
    
clear Foi Toi Tbase
 
 
Foi=[20 28];
Toi=[0.2 0.35];
Tbase=[-0.3 -0.15];
    
    
    
    
for n_sub = 12:length(subject)
    
    %% COV_noise

    infile1 = fullfile(raw_exp,'raw',[subject{n_sub} '_raw_tsss.fif']);
 
    oridata_obj=
    load('S13_social.mat');  % pre=obj

    
    cfg        = [];
    cfg.channel      = 'MEG';
    data1   = ft_selectdata(cfg, HS);
    
    cfg        = [];
    cfg.channel      = 'MEG';
    data2   = ft_selectdata(cfg, LS);
    
    
      %% baseline for obj
     
     COV_noise = [];
     COV_noise.COV = zeros(306);     
     keep_noise = find(data1.time{1} > min(Tbase) & data1.time{1} < max(Tbase));  % baseline boundaries defined here
     for n_trial = 1:length(data1.trial)
            COV_noise.COV = COV_noise.COV + data1.trial{n_trial}(:,keep_noise)*data1.trial{n_trial}(:,keep_noise)';
     end
     COV_noise.COV = COV_noise.COV/length(keep_noise)/n_trial;
    
    %%  baseline for verb
     COV_noise2 = [];
     COV_noise2.COV = zeros(306);     
     keep_noise2 = find(data2.time{1} > min(Tbase) & data2.time{1} < max(Tbase));  % baseline boundaries defined here
     for n_trial = 1:length(data2.trial)
            COV_noise2.COV = COV_noise2.COV + data2.trial{n_trial}(:,keep_noise2)*data2.trial{n_trial}(:,keep_noise)';
     end
     COV_noise2.COV = COV_noise2.COV/length(keep_noise2)/n_trial;
     
   
     fwdfile = fullfile(subjects_dir,subject{n_sub},'meg',[subject{n_sub},'_from_MNI-Mirror-5-src-fwd.fif'])
     fwd = CM_lf2comp(fwdfile,COV_noise)

%% CSD at the time of interest
    % CSD at the time of interest  
  
         cfg        = [];
    cfg.channel      = 'MEG';
    cfg.latency = Toi;
    data_obj   = ft_selectdata(cfg, HS);
    data_verb  = ft_selectdata(cfg, LS);
    cfg        = [];
    cfg.channel      = 'MEG';
    cfg.latency = Tbase;
    data_bsl_obj   = ft_selectdata(cfg, HS);
    data_bsl_verb  = ft_selectdata(cfg, LS);
   
    cfg = [];
    cfg.channel   = 'MEG';
    cfg.method    = 'mtmconvol';
    cfg.output    = 'powandcsd';
    cfg.taper =  'hanning';
    cfg.foi  = 20:1:21  ;
    cfg.toi = 0.2:0.01:0.35;
    cfg.t_ftimwin =ones(length(cfg.foi),1).*0.5;
    freq_source_obj = ft_freqanalysis(cfg,  HS);
    freq_source_obj = CM_reshape_ft_CSD(freq_source_obj);
    X = find(isnan(freq_source_obj.CSD));
    cfg = [];
    cfg.method    = 'mtmconvol';
    cfg.channel   = 'MEG';
    cfg.output    = 'powandcsd';
    cfg.taper =  'hanning';
    cfg.foi  = 20:1:21  ;
    cfg.toi = -0.35:0.01:-0.2;
    cfg.t_ftimwin = ones(length(cfg.foi),1).*0.5;
    freq_source_base_obj = ft_freqanalysis(cfg, HS);
    freq_source_base_obj = CM_reshape_ft_CSD(freq_source_base_obj);
    X = find(isnan(freq_source_base_obj.CSD));
    
    
   %% data_verb
    cfg = [];
    cfg.channel   = 'MEG';
    cfg.method    = 'mtmconvol';
    cfg.output    = 'powandcsd';
    cfg.taper =  'hanning';
    cfg.foi  = 20:1:21  ;
    cfg.toi = 0.2:0.01:0.35;
    cfg.t_ftimwin = ones(length(cfg.foi),1).*0.5;
    freq_source_verb = ft_freqanalysis(cfg, LS);
    freq_source_verb = CM_reshape_ft_CSD(freq_source_verb);
    X = find(isnan(freq_source_verb.CSD));
    
    cfg = [];
    cfg.channel   = 'MEG';
    cfg.method    = 'mtmconvol';
    cfg.output    = 'powandcsd';
    cfg.taper =  'hanning';
    cfg.foi  = 20:1:21  ;
    cfg.toi = -0.35:0.01:-0.2;
    cfg.t_ftimwin = ones(length(cfg.foi),1).*0.5;
    freq_source_base_verb = ft_freqanalysis(cfg, LS);
    freq_source_base_verb = CM_reshape_ft_CSD(freq_source_base_verb);
    X = find(isnan(freq_source_base_verb.CSD));
  
    
    CSD_pooled = [];
    CSD_pooled.infile = infile1;
    CSD_pooled.bads = [];
    CSD_pooled.CSD = real(mean(freq_source_obj.CSD(:,:,:),3) + mean(freq_source_base_obj.CSD(:,:,:),3))/2;
    Tinv.fwdfile = fwdfile;
    Tinv = CM_sf_LCMV(Tinv,CSD_pooled);
   
    CSD_pooled = [];
    CSD_pooled.infile = infile1;
    CSD_pooled.bads = [];
    CSD_pooled.CSD = real(mean(freq_source_verb.CSD(:,:,:),3) + mean(freq_source_base_verb.CSD(:,:,:),3))/2;
    Tinv2.fwdfile = fwdfile;
    Tinv2 = CM_sf_LCMV(Tinv2,CSD_pooled);
    
    
    
%% source power pre
% Estimate the power in the source space (left and right multiplication
% of the CSD by the inverse operator) first for the time of interest
% and then for the baseline.
      Nvox = double(Tinv.nsource);
    % for the time of interest
      P = zeros(1,Nvox);
      T_CSD = Tinv.Tinv*mean(freq_source_obj.CSD(:,:,:),3);
      for n = 1:Nvox
          Cxx = T_CSD([-1 0]+2*n,:)*Tinv.Tinv([-1 0]+2*n,:)';
          P(n) = max(abs(eig(Cxx)));
      end
    % for the baseline
      Pbase = zeros(1,Nvox);
      T_CSD = Tinv.Tinv*mean(freq_source_base_obj.CSD(:,:,:),3);
        for n = 1:Nvox
            Cxx = T_CSD([-1 0]+2*n,:)*Tinv.Tinv([-1 0]+2*n,:)';
            Pbase(n) = max(abs(eig(Cxx)));
        end
        
    obj = P
    obj_base = Pbase
%%

%% source power post
% Estimate the power in the source space (left and right multiplication
% of the CSD by the inverse operator) first for the time of interest
% and then for the baseline.
      Nvox = double(Tinv2.nsource);
    % for the time of interest
      P = zeros(1,Nvox);
      T_CSD = Tinv2.Tinv*mean(freq_source_verb.CSD(:,:,:),3);
      for n = 1:Nvox
          Cxx = T_CSD([-1 0]+2*n,:)*Tinv2.Tinv([-1 0]+2*n,:)';
          P(n) = max(abs(eig(Cxx)));
      end
    % for the baseline
      Pbase = zeros(1,Nvox);
      T_CSD = Tinv2.Tinv*mean(freq_source_base_verb.CSD(:,:,:),3);
        for n = 1:Nvox
            Cxx = T_CSD([-1 0]+2*n,:)*Tinv2.Tinv([-1 0]+2*n,:)';
            Pbase(n) = max(abs(eig(Cxx)));
        end
    
      verb = P
      verb_base = Pbase
%% Estimate the power ratio and take care of that some sources are ouside 

   
      fwd = mne_read_forward_solution(fwdfile)
%     
      Pratio_verb =zeros(1,fwd.src.np);   % change the fwd
      Pratio_verb(fwd.src.vertno) = (verb - verb_base)./verb_base   ;
      
      Pratio_obj =zeros(1,fwd.src.np);   % change the fwd
      Pratio_obj(fwd.src.vertno) = (obj - obj_base)./obj_base ;
     

    % convert the vector Pratio into a 3D overlay corregistered to MNI MRI
      map_mri_verb = CM_overlay_MNI((Pratio_verb)',5,subjects_dir);   
      map_mri_verb(map_mri_verb==0) = 0; % set the power ratio to 1 for voxels outside the BEM
%     
      map_mri_obj = CM_overlay_MNI((Pratio_obj)',5,subjects_dir);   
      map_mri_obj(map_mri_obj==0) = 0; % set the power ratio to 1 for voxels outs
%       
    %% write down the overlay. Overlays starticlearng with 'w' are normalized to
    % MNI template. Those starting with 'iw' are normalized to subject's
    % anatomy and should be visualized on the 'orig' MRI available in the
    S = spm_vol(fullfile(subjects_dir , 'MNI' , 'norm' , 'MNI.nii'));
    if not(exist(fullfile(subjects_dir , subject{n_sub} , exp)))
        mkdir(fullfile(subjects_dir , subject{n_sub} , exp));
    end
    

    outputfile = 'social'; 
    map_mri = map_mri_verb;
    
    S.fname = fullfile(subjects_dir, outputfile,[subject{n_sub} '_LS_beta.nii']);
    S.dt(1) = 4;
    S.pinfo(1) = max(map_mri(:))/(2^15-1);
    S.pinfo(2) = 0;
    S.pinfo(3) = 0;
    S.private.dat.fname = S.fname;
    S = spm_write_vol(S,map_mri);

    
    clear map_mri
    
 
    map_mri = map_mri_obj;
    
    S.fname = fullfile(subjects_dir, outputfile,[subject{n_sub} '_HS_beta.nii']);
    S.dt(1) = 4;
    S.pinfo(1) = max(map_mri(:))/(2^15-1);
    S.pinfo(2) = 0;
    S.pinfo(3) = 0;
    S.private.dat.fname = S.fname;
    S = spm_write_vol(S,map_mri);
        
    clear map_mri
       
 end
 
 
