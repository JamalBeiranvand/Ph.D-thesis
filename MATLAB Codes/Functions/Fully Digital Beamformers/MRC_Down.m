 % Author  : Jamal Beiranvand (Jamalbeiranvand@gmail.com)
 % google scholar: https://scholar.google.com/citations?user=S6LywwsAAAAJ&hl=en
 % Date : 29/05/2020
function Fmrc=MRC_Down(H)
     A=H';
     alpha=1/trace(A*A');
     Fmrc=A*sqrt(alpha);
end