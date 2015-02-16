
function G=geometric_moments_orig(X,K,N,num_facets,num_vertices)
% Computes the geometric moments
% of the volumetric object given

% PRE-COMPUTATIONS
display('Trinomial')
% Matrix to store the pre-computed trinomial
tri_matrix=trinomial_matrix(N);
% Volume of each tetrahedra
Vol=zeros(num_facets,1);
% Matrix that stored the monomial combinations
%C_temp=zeros(num_vertices,N+1,N+1,N+1);
%keyboard
%size(C)

display('Dabc_orig')
%matlabpool close
%matlabpool local 8
D=zeros(num_facets,N+1,N+1,N+1);
C1=zeros(num_facets,N+1,N+1,N+1);
for facet=1:num_facets
    % pre-compute D
    tic
    if mod(facet,1000) == 0
        display(['Facet: ' num2str(facet) ' of ' num2str(num_facets)])
    end
    [C_temp Vol_temp] = D_CV_orig(facet,num_vertices,N,X,K,tri_matrix);
    D(facet,:,:,:)=Dabc_orig(C_temp,N);
    C1(facet,:,:,:) = C_temp(1,:,:,:);
    Vol(facet,1) = Vol_temp;
    if mod(facet,1000) == 0
        display(['Computed facet: ' num2str(facet) ' in ' num2str(toc) 'seconds'])
    end
end
clear C_temp Vol_temp

display('D_SG_orig')
%matlabpool close
%matlabpool local 8
% GEOMETRIC MOMENTS;
F = factorial_precalc(N);
G=zeros(N+1,N+1,N+1);
for i=0:N
    tic
    display(['Order: ' num2str(i) ' of ' num2str(N)])
    G(i+1,:,:)=D_SG_orig(num_facets,i,N,C1,D,Vol,F);
    display(['Computed Order ' num2str(i) ' in ' num2str(toc) 'seconds'])
end







