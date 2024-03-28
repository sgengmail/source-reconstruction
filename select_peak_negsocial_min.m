exp_dir = '/bcbl/home/home_n-z/sgeng/MEGtool/toolboxes/';

% tbpath = '~/Dropbox/toolboxes/'
toolboxes = {'cartographie_motrice' 'fieldtrip-20170911'};
for tool = toolboxes
    addpath(genpath(fullfile(exp_dir,tool{1})))  
end  
    
  subjects = {
        'S1'
        'S2';
        'S3';
         'S4'
        'S5';
        'S6';
        'S7';
        'S8';
        'S9';
         'S11'
        'S12';
        'S13';
        }; 
  
raw_exp = 'G:\Premeg\shuang\Freesurfer';   
meg_dir = 'G:\Premeg\shuang\MEG';
cd(fullfile(raw_exp));
subjects_dir = fullfile(raw_exp);
raw_dir =  fullfile(raw_exp,'raw');
exp = 'social';
group_fold = fullfile(subjects_dir,exp)


Nave = length(subjects);
spm_file = [];
for k = 1:Nave
    spm_file{k} = fullfile(group_fold,[subjects{k} '_HS_beta.nii']);
end

Nave = length(subjects);
spm_file_perm_stat = [];
for k = 1:Nave
    spm_file_perm_stat{k} = fullfile(group_fold,[subjects{k} '_LS_beta.nii']);
end

V = spm_vol(spm_file{1});
Z = zeros([prod(V.dim) Nave]);
for k = 1:Nave
    V = spm_vol(spm_file{k});
    C = log(spm_read_vols(V)+1);
    Vp = spm_vol(spm_file_perm_stat{k});
    Cp = log(spm_read_vols(Vp)+1);
    Z(:,k) = (C(:)-Cp(:))/2;
end


NaNs = sum(isnan(Z),2);
one_NaN = find(NaNs == 1);
subZ = Z(one_NaN,:);
subZ(isnan(subZ)) = 0;
Z(one_NaN,:) = subZ;

% sub_ind to speed up the processing
sub_ind = ones(V.dim);
set_to_0 = 1:V.dim(1);
set_to_0(3:5:end) = [];
sub_ind(set_to_0,:,:) = 0;
sub_ind(:,set_to_0,:) = 0;
sub_ind(:,:,set_to_0) = 0;
sub_ind = find(sub_ind);
to_rm = find(max(abs(Z(sub_ind,:)),[],2) == 0);
sub_ind(to_rm) = [];

% choose the permutation and start the loop for the permutation test
Nsim = 1000;
rand('seed',sum(clock)*1000)
stat = zeros(1,Nsim+1);
for k = 1:Nsim+1
    if ~mod(k,20)
        disp(['Processing simulation ' num2str(k) '/' num2str(Nsim)])
        pause(0.0001)
    end
   
    % determine permutations
    perm = sign(randn(Nave,1));
    if k == (Nsim + 1)
        perm = ones(Nave,1);
    end
   
    % Extract the statistic
    sumZ = Z(sub_ind,:)*perm;
    stat(k) = min(sumZ/Nave);    %%%%enlever abs si mesure l'un moins l'autre%%%
end

Zave = Z*perm/Nave;
p = mean(min(Zave) > stat(1:Nsim));
Th = prctile(stat(1:Nsim),01);
Zave = reshape(Zave,V.dim);


% write zmap. this part needs to be re-written
W = V;
W.dt(1) = 4;
W.pinfo(1) = (max(Zave(:))-min(Zave(:)))/(2^15-1);
W.pinfo(2) = min(Zave(:));
W.pinfo(3) = 0;
if ~exist(fullfile(group_fold,'group_ave'))
    unix(['mkdir ' fullfile(group_fold,'group_ave')])
end
W.fname = fullfile(group_fold,'group_ave',['min' '_Th_' num2str(Th) '.nii']);
W.private.dat.fname = W.fname;
W = spm_write_vol(W,Zave);
% 
%% transform form Z to T



%write a p-valed map
sort_stat = sort(stat(1:Nsim),'descend');
pmap_one = ones(size(Zave));
keep = find(Zave > Th);
pmap_one(keep) = 0;
for ind = 1:Nsim*0.05;
    pmap_one(keep) = pmap_one(keep) + (abs(Zave(keep)) >= sort_stat(ind))/Nsim;
end
pmap = pmap_one;
pmap(pmap == 0) = 0.001;
pmap = -log10(pmap);


Va = W;
Va.fname = fullfile(group_fold,'group_ave',['ventral''ventral_pmap.nii']);
Va.private.dat.fname = Va.fname;
Va.pinfo(1) = 3/(2^15-1);
Va.pinfo(2) = 0;
Va.pinfo(3) = 0;
Va = spm_write_vol(Va,pmap);


% smooth the map
FWHM = 8;
kern = CM_gaussian_kern(FWHM,W.dim);
kern = kern/sum(kern(:));
Fkern = fftn(kern);
FZave = fftn(Zave);
FZave = FZave.*Fkern;
Zave = real(ifftn(FZave));

[cluster_neg,d] = CM_local_maxima(W,abs(Zave),0.03);  % threshold of 0.13 fixed based on inspection of the map (to uncover the left SM1 source)

% manual selection
keep = [1:10]; % manual selection (3 first are significant, last is added)
cluster_neg = cluster_neg(keep);
matfile = fullfile(group_fold,'group_ave',['cluster_group_pred']);
save(matfile,'cluster_neg')

clear rr_group_MNI

% store the coordinates

cluster=cluster_neg;

coord_group = [];
rr_coord_group = [];
rr_group = [];
for n_source = 1:length(cluster);
    coord_group(:,n_source) = mean(cluster(n_source).mni,1);
    rr_group_MNI(:,n_source) = round(mean(cluster(n_source).mni,1));
    rr_group(:,n_source) = round(mean(cluster(n_source).vox,1));
end

mat_free2spm = eye(4); mat_free2spm(1:3,4) = [0 -16 20]; % see CM_overlay for the computation

vals = zeros(length(cluster), length(subjects), 2);
pvals = zeros(length(cluster),1);
stats = cell(length(cluster),1);
for n_source = 1:length(cluster);
    for n_sub = 1:length(subjects)
        to = subjects{n_sub};
        from = 'MNI';
        file1 = fullfile(group_fold,[subjects{n_sub} '_HS_beta.nii']);
        file2 = fullfile(group_fold,[subjects{n_sub} '_LS_beta.nii']);
       
        V1 = spm_vol(file1);
        V2 = spm_vol(file2);
        C1 = spm_read_vols(V1);
        C2 = spm_read_vols(V2);
        Z = (C1(:)+C2(:))/2;
        Nii_mat = reshape(Z,V1.dim);
        clear Z;
       
        % smooth the maps
        FWHM = 5;
        kern = CM_gaussian_kern(FWHM,V1.dim);
        kern = kern/sum(kern(:));
        Fkern = fftn(kern);
       
        FNii_mat = fftn(Nii_mat);
        FNii_mat = FNii_mat.*Fkern;
        Nii_mat = real(ifftn(FNii_mat));
        Z = Nii_mat;
       
        C1_mat = fftn(C1);
        C1_mat = C1_mat.*Fkern;
        D1 = real(ifftn(C1_mat));
       
        C2_mat = fftn(C2);
        C2_mat = C2_mat.*Fkern;
        D2 = real(ifftn(C2_mat));
       
        coord = [];
        %for n_source = 1:length(cluster);
        tmp = rr_group(:,n_source);   %% coordinates of 10 sourcess 
        r = 7;
        d = 2*r+1;
        x = tmp(1)-r:tmp(1)+r; x = repmat(x',[1 d d]);                 x = x(:);
        y = tmp(2)-r:tmp(2)+r; y = repmat(y,[d 1 d]);                  y = y(:);
        z = tmp(3)-r:tmp(3)+r; z = repmat(reshape(z,[1 1 d]),[d d 1]); z = z(:);  %% XYZ
        
        rr = [x y z];
        keep = find(sum((rr - repmat(tmp',[d^3 1])).^2,2)<r^2);
        rr = rr(keep,:);
        n = CM_xyz2n(rr,size(Z));  %Convert (x y z) coordinates of a 3D matrix of size sizemat into the    index of this matrix corresponding to (x y z)
        %% n shi di n ge voxel
        
       [M,P] = max(Z(n));   % P index of maximum intensity, M the value of intensity  
        E1 = D1(n);
        E2 = D2(n);
        %coord(:,n_source) = inv(mat_free2spm)*V1.mat*[rr(P,:) 1]';
        vals(n_source, n_sub, :) = [E1(P) E2(P)];
    end
    %log_vals = log(vals);
   
    [h,pvals(n_source),ci,stats{n_source}] = ttest(vals(n_source,:,1), vals(n_source,:,2));
end
