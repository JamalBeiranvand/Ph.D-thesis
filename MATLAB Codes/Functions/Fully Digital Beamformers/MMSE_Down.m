 % Author  : Jamal Beiranvand (Jamalbeiranvand@gmail.com)
 % google scholar: https://scholar.google.com/citations?user=S6LywwsAAAAJ&hl=en
 % Date : 29/05/2020
 function Fmmse=MMSE_Down(H,SNR)
 [K,~]=size(H);
 A=H'/(H*H'+K/SNR*eye(K));
 alpha=1/trace(A*A');
 Fmmse=A*sqrt(alpha);
 end