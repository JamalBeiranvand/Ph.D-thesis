 % Author  : Jamal Beiranvand (Jamalbeiranvand@gmail.com)
 % google scholar: https://scholar.google.com/citations?user=S6LywwsAAAAJ&hl=en
 % Date : 29/05/2020
function Fzf=ZF_Down(H)
         A=H'/(H*H');
         alpha=1/trace(A*A');
         Fzf=A*sqrt(alpha);
end
