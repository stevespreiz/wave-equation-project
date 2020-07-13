function [unp1] = timeStep(sigma,ja,jb,unm1,un,unp1,nD,oacc)
% Main time step over middle values

if nD == 1
    for i = ja:jb
        unp1(i) =  2*un(i) - unm1(i) + sigma^2*(un(i+1)-2*un(i)+un(i-1));
        if oacc > 2
            unp1(i) = unp1(i) - (sigma^2-sigma^4)/12*(un(i+2)-4*un(i+1)+6*un(i)-4*un(i-1)+un(i-2));
            if oacc > 4
                unp1(i) = unp1(i) + (sigma^2/90-sigma^4/72+sigma^6/720)* ...
                    (un(i+3)-6*un(i+2)+15*un(i+1)-20*un(i)+15*un(i-1)-6*un(i-2)+un(i-3));
%                                 unp1(i) = 2*un(i) - unm1(i)+...
%                     sigma^2/180*(2*un(i-3)-27*un(i-2)+270*un(i-1)-490*un(i)+270*un(i+1)-27*un(i+2)+2*un(i+3))+...
%                     sigma^4/12/6*(-1*un(i-3)+12*un(i-2)-39*un(i-1)+56*un(i)-39*un(i+1)+12*un(i+2)-1*un(i+3))+...
%                     sigma^6/720*(un(i+3)-6*un(i+2)+15*un(i+1)-20*un(i)+15*un(i-1)-6*un(i-2)+un(i-3));
            end
        end
    end
end

end

