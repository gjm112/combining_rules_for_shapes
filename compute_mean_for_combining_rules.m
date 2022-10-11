
cd /Users/gregorymatthews/Dropbox/full_shape_classification_SRVF/code/matlab/

%Read in data
teeth = readtable("/Users/gregorymatthews/Dropbox/combining_rules_for_shapes/teeth_out.csv")

%teethBWgladysvale500matrix = renamevars(teethBWgladysvale500matrix,["Var1"],["image"])


%get the number of rows and cols
n_rows = size(teeth,1);
n_cols = size(teeth,2); 

%blank array storage for matric for each image
teeth_data = zeros(2,n_cols,n_rows/2);

%loop to assign each tooth's coordinates to the appropriate spot in the
%array
for i=1:2:n_rows
    %get index for the image we are on 
    j = find((1:2:n_rows)==i);
    
    %assign matrix for image j 
     X = teeth{i:(i+1), 2:n_cols};
    
    %resample so all of the curves have 100 points
    teeth_data(:,:,j) = ReSampleCurve(X,n_cols);
end


%Compute mean
[out_q,q,E] = FindElasticMeanFast(teeth_data)
out_beta = q_to_curve(out_q)



%return it to R
csvwrite("/Users/gregorymatthews/Dropbox/combining_rules_for_shapes/mean_q.csv", out_q)
csvwrite("/Users/gregorymatthews/Dropbox/combining_rules_for_shapes/mean_beta.csv", out_beta)


cd /Users/gregorymatthews/Dropbox/combining_rules_for_shapes/
save('out_q')
save('out_beta')


%now project teeth to tanget space of the mean.
%numPCs = 10
%[VV PC] = FindTangentFeatures(out_q,teeth_data,numPCs)

%csvwrite("/Users/gregorymatthews/Dropbox/combining_rules_for_shapes/teeth_pop_coords_in_tan_space_of_overall_mean.csv",VV)
%csvwrite("/Users/gregorymatthews/Dropbox/combining_rules_for_shapes/teeth_pop_PC_in_tan_space_of_overall_mean.csv", PC)






    
