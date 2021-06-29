function U = laplacian_surface_editing_2D(vertex,BI,BC)

%The file is an implementation of 'Laplacian surface editing' in 2D. 
%Inputs: vertex 
  %   vertex  #vertex by 2 list of rest domain positions
  %   bi       #b list of indices of constraint (boundary) vertices
  %   bc      #b by 2 list of constraint positions for b
  
%Output: 
  %   U       #V by dim list of new positions

% By raymond @ smartee on 28/06/2021

n = length(vertex);

% the Laplacian matrix (uniform weighting)
L = spdiags(ones(n,1),0,n,n) - spdiags(ones(n,1),1,n,n);
L = L+L';
L(1,n)= -1;
L(n,1) = -1;
L = L./2;


delta = L*vertex;


% we want to construct the matrix of the system for v-primes
L_prime = [   L     zeros(n)   % the x-part
	       zeros(n)    L    ]; % the y-part


 for i = 1:n
     ring = [(1 + mod(i-2,n)), i, (1 + mod(i,n))];
      V = vertex(ring,:)';
      V = [V
      ones(1,length(ring))];%Here is ones matrix, multiplying V' becomes A in formula (10). Such writing is for associating multiplication factors of v'. 
   	  %The first row of V is x part,the second row is y part, the third one is z part, the elements in last row are all ones.
        C = zeros((length(ring)-1) * 3, 4);
   % ... Fill C in
  for r=1:length(ring)
    C(r,:) =                [V(1,r)       V(2,r)  V(3,r)      0  ];
    C(length(ring)+r,:) =   [V(2,r)  (-1)*V(1,r)       0  V(3,r) ];
  end;  
   Cinv = pinv(C);
  s =  Cinv(1,:);
  a =  Cinv(2,:);

 
  delta_i = delta(i,:)';
  delta_ix = delta_i(1);
  delta_iy = delta_i(2);  
  
   % T*delta gives us an array of coefficients  
   % T*delta*V' equals to T(V')*delta in formula (5)
  Tdelta = [delta_ix*s      + delta_iy*a 
	        delta_ix*(-1)*a + delta_iy*s];
        
  % updating the weights in Lx_prime, Ly_prime
  % Note that L_prime has already containted L. Here L_prime represents T(V')*delta - L(V') in formula(5)
  L_prime(i,[ring (ring + n)]) = L_prime(i,[ring (ring + n)]) +...
                                              (-1)*Tdelta(1,:);
  L_prime(i+n,[ring (ring + n)]) = L_prime(i+n,[ring (ring + n)]) +...
                                                (-1)*Tdelta(2,:);
 end
   
% weight for the constraints
w=1;

% building the least-squares system matrix
A_prime = L_prime;
rhs = zeros(2*n,1);




for j=1:length(BI)
  A_prime = [A_prime
	     w*((1:(2*n))==BI(j))
	     w*((1:(2*n))==(BI(j)+n))];
         
  rhs = [rhs
	 w*BC(j,1)
	 w*BC(j,2)];
end;

% solving for v-primes
xyz_col = A_prime\rhs;
U = [xyz_col(1:n) xyz_col((n+1):(2*n))];