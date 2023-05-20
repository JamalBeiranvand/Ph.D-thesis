% Author  : Jamal Beiranvand (Jamalbeiranvand@gmail.com)
 % google scholar: https://scholar.google.com/citations?user=S6LywwsAAAAJ&hl=en
 % Date : 29/05/2020
 function Wmmse=MMSE_Up(H,SNR)
         [~,K]=size(H);
         Wmmse=H/(H'*H+1/SNR*eye(K));
 end