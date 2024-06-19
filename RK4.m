function [x,y]=RK4(fid,f,a,eta,h,n,sol)
%funcion que implementa el metodo de runge kutta clásico
y=eta;
x=a;
maxerr=0;
error=0;
for i=1:n
    if (nargin==7)
        %maxerr=0;
        %error=0;
        escribe_cabecera(fid,x,y,error);
        escribe_paso(fid,i,x,y,error);
        error=y-sol(x);
        if norm(y-sol(x))>maxerr
            maxerr=norm(y-sol(x));
        end
    else
        escribe_cabecera(fid,x,y);
        escribe_paso(fid,i,x,y);

    end 
    %%%%terminos intermedios%%%%%
    %xn1=x;
    xn2=x+h/2;
    xn3=xn2;
    xn4=x+h;
    kn1=f(x,y);
    kn2=f(xn2,y+kn1*(h/2));
    kn3=f(xn3,y+kn2*(h/2));
    kn4=f(xn4,y+kn3*h);
    %%%%actualizacion terminos%%%%    
    y=y+(h/6)*(kn1+2*kn2+2*kn3+kn4);
    x=x+h;
end
