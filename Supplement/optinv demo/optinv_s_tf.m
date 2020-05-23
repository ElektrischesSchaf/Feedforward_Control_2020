% 	optimal inversion function
%	Usage: [uopt,yopt] = optinv01(num,den,r,q,omega,y,t)
%	
%	Notes
%
%	num(s)/den(s) describes a stable transfer function
%	system must be non-hyperbolic, i.e., no zeros on the imaginary axis
%	r and q are defined as a function of omega, i.e. specify frequencies 
%	values in omega and the values of r and q for those frequencies
%	If r and q are not specified for a frequency then r =0 and q=1; 
%	t row vector, must be evenly spaced time points (and has an odd length)
%	it is advisable that y be zero at beginning and end.
%	y must be same size as time vector "t"
%	the output uopt is the optimal inverse-input into the system
%	yopt is the modified output
function [uopt,yopt] = optinv_s_tf(num,den,r,q,omega,y,t)

%Related Papers 
%[1] 	Q. Zou and S. Devasia "Preview-based Stable-Inversion for Output Tracking," 
%	ASME J. of Dynamic Systems, Measurement and Control, 
%	Vol. 121 (4), pp. 625-630, December 1999. 
%[2]	S. Devasia "Approximated Stable Inversion for Nonlinear Systems with 
%	Nonhyperbolic Internal Dynamics," IEEE Trans. on Automatic Control, 
%	Vol. 44 (7), pp. 1419-1425, July 99.
%[3]	J.S. Dewey, K. K. Leang and S. Devasia "Experimental and Theoretical Results 
%	in Output-Trajectory Redesign for Flexible Structures," ASME Journal of 
%   	Dynamic Systems, Measurement, and Control, Vol. 120 (4), pp. 456-461, Dec. 1998. 
%[4]  	S. Devasia and B. Paden "Stable Inversion for Nonlinear Nonminimum-Phase 
%	Time-Varying Systems," IEEE Trans. on Automatic Control, 
%     	Vol. 43 (2), pp. 283-288, Feb. 1998. 
%[5]	S. Devasia, B. Paden and C. Rossi "Exact-Output Tracking Theory for Systems 
%	with Parameter Jumps," International Journal of Control, 
%	Vol. 67 (1), pp. 117-131, May 1997.       
%[6]	S. Devasia, D. Chen and B. Paden "Nonlinear Inversion-Based Output Tracking,"
%	IEEE Transactions on Automatic Control, Vol. 41 (7), pp. 930-942, July 1996.    
%[7] 	D. Croft and S. Devasia "Vibration Compensation for High Speed Scanning 
%	Tunneling Microscopy," Review of Scientific Instruments published by 
%	the American Institute of Physics, Vol. 70 (12), pp. 4600-4605, December 1999.
  
delt 		= t(2)-t(1);
Nt 		    = length(t);
Nf 		    = (Nt-1)/2;
del_freq    = ((2*pi/delt)/2)/Nf;
i 			= sqrt(-1);
u 			= zeros(size(t));
yopt 		= zeros(size(t));
yf 		    = fft(y);
ome_max 	= max(omega); 
%ym = abs(yf); yp=(180/pi)*unwrap(angle((yf))); 
%ww=(1:1:Nt)*(2*pi/delt);plot(ww,abs(ym))

for jj=1:1:Nf,
   ome	= del_freq*jj; 
   w  	= i*ome;
   if ome < ome_max
      rr = interp1(omega,r,ome); 
      qq = interp1(omega,q,ome); 
   else
      rr = 1; 
      qq = 0;
   end
   nn 			= polyval(num,w)*polyval(num,-w); 
   dd 			= polyval(den,w)*polyval(den,-w);
   dw 			= polyval(den,w); nw = polyval(num,w);
   yopt(jj+1) 	= (nn*qq/(nn*qq + dd*rr))*yf(jj+1);
   u(jj+1) 		= (dw/nw)*yopt(jj+1);
   jjj 			= Nt+1-jj;
   yopt(jjj) 	= conj(yopt(jj+1));
   u(jjj) 		= conj(u(jj+1));
end

% the DC gain effects
rr 		= interp1(omega,r,0); 
qq 		= interp1(omega,q,0); 
nn 		= polyval(num,0)*polyval(num,-0);
dd 		= polyval(den,0)*polyval(den,-0);
dw 		= polyval(den,0);
nw 		= polyval(num,0);
yopt(1) 	= (nn*qq/(nn*qq + dd*rr))*yf(1);
u(1) 		= (dw/nw)*yopt(1);

yopt		= real(ifft(yopt));
uopt		= real(ifft(u));

return