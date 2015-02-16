%add the path of the zernike code
addpath('../matlab_only'); rmpath('../mex')
%addpath('../mex'); rmpath('../matlab_only')

[V,F] = read_vtk('Parallelepiped.vtk'); V = V'; F = F';

num_vertices = 3;
num_facets   = size(F,1);
N = 20; %code will only work upto order 20 (can be fixed later)

% ZERNIKE MOMENTS
G=geometric_moments_orig(V,F,N,num_facets,num_vertices);
Z=zernike(G,N);
Descriptors=feature_extraction(Z,N);
ZM= reformat_zernike(Z,N);

%compare to your code
load ZMvtk
max_abs_diff = max(abs(abs(ZM) - abs(ZMvtk))); %error should be around 1e-8
display(['Max Absolute Difference: ' num2str(max_abs_diff) ]) 