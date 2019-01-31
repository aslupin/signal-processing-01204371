t = 10;
a1 = 13;
a2 = 20;
converterToRadiasn = pi/180;
phase1 = 88 * converterToRadiasn;
phase2 = -30 * converterToRadiasn;
x1 = a1 * cos(2*pi*(3000 * t) + phase1)
x2 = a2 * cos(2*pi*(3000 * t) + phase2)
plot(x1)
plot(x2)

% function	Z = replacez(A)
% %REPLACEZ	Function that replaces the negative elements
% %	of a matrix with the number 77
% %	usage:
% %	Z = replacez(A)
% %	A  = input matrix whose negative elements are to
% %	be replaced with 77
% %
% [M,N] = size(A);
%     i = 1;
%     % for i=1:M
%     while(i <= M)
%         j = 1;
%         % for j=1:N
%         while(j <= N)
%             Z(i,j) =  (A(i,j) < 0) * 77 + (A(i,j) >= 0) * A(i,j);
            
%             % if A(i,j) < 0
%             %     Z(i,j) = 77;
%             % end
%             % else
%             %     Z(i,j) = A(i,j);
%             % end
%             j = j + 1;
%         end 
%         i = i + 1;
%     end
% end



% A = randn(6,3)
% % สร้าง matrix ขนาด 6 * 3 โดย แรนดอมค่าแบบ normal distribution
% % B = (A>0)
% A = A .* (A>0) 
% % .* คือ operation ในการคูณ array โดย operand ด้านขวาคือการ map array A เป็น 1 หรือ 0 โดนมี condition ว่า A มากกว่า 0 หรือไม่


% function Z = expand(xx,ncol)
% %EXPAND	Function to generate a matrix Z with identical columns
% %	equal to an input vector xx
% %	usage:
% %	Z = expand(xx,ncol)
% %	xx = the input vector containing one column for Z
% %	ncol = the number of desired columns
% %
%     % xx = xx(:);	%-- makes the input vector x into a column vector
%     % Z = zeros(length(xx),ncol);
%     % % for i=1:ncol
%     % Z(:,1:1:ncol) = xx;
%     xx = xx(:);	%-- makes the input vector x into a column vector
%     Z = zeros(length(xx),ncol);
%     i = 1
%     while(i <= ncol)
%     % for i = 1:ncol
%         Z(:,i) = xx;
%         i =  i + 1;
%     end
% end

% yy = ones(1,10)
% a = expand(yy,4)


% tt = -2 :0.05 :	3;
% xx = sin( 2*pi*0.789*tt );
% plot( tt, xx ), grid on %<--- plot a sinusoid title(’TEST PLOT of SINUSOID’)
% plot(0.5*cos( 2*pi*0.789*tt ))


% SOUND
% dur = 1.0;
% fs = 8000;
% tt = 0: (1/fs):dur; 
% xx = sin(2*pi*2000*tt ); 
% length(xx)
% sound( xx, fs )
% soundsc(xx, fs)

% Function

% function xx = cosgent(f,dur)
    
    
% 	tt = [0:1/(20*f):dur];
% 	yy = cos(2*pi*f*tt);
%     xx = yy;
%     return;
% end
% a = cosgent(2,20)


% function y = average(x)
% if ~isvector(x)
%     error('Input must be a vector')
% end
%     y = sum(x)/length(x); 
% end
% z = 1:99;
% average(z)

% function xx = cosgen(f,dur)
% 	%COSGEN	Function to generate a cosine wave
% 	%	usage:
% 	%	xx = cosgen(f,dur)
% 	%	f = desired frequency
% 	%	dur = duration of the waveform in seconds
% %
% 	tt = [0:1/(20*f):dur];	%  gives 20 samples per period 
% 	yy = cos(2*pi*f*tt);
%     xx = yy;
% end


% function [sum,prod] = sumprod(x1,x2)
% %SUMPROD	Function to add and multiply two complex numbers
% %	usage:
% %	[sum,prod] = sumprod(x1,x2)
% %	x1 = a complex number
% %	x2 = another complex number
% %	sum  = sum of x1 and x2
% %	prod = product of x1 and x2
% %

% sum = x1+x2;
% prod = x1*x2;
% end
% [a,b] = sumprod(1,1)

% yy = ones(7,1) * rand(1,4) % [ 7 * 1] x [ 1 * 4] 
% xx = randn(5,5)
% jj = ones(6,1)
% xx(1,:)

% yy = xx(ones(6,1),:)


