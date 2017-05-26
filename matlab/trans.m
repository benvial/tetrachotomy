function out=trans(Z,alpha0,norm)
c=299792458;

global op or einf

% thickness of the slab (can be a multilayer with N layers then h is a vector of size N)
h=[3]*1e-6;

pola = 'TE';
Omega=[];
for ilr=1:size(Z,1)
    for ili=1:size(Z,2)
        
        clear C s d T M
        Omega=Z(ilr,ili);
        lambda0=2*pi*c./Omega;
        
        % dispersive epsilon
        ep=einf-op^2./(Omega.^2-or^2);
        
        
        % eps contains the permitivities of the different layers
        % [superstrate, layer 1, layer 2,... layer N, substrate]
        eps1=1; 
        eps2=ep;
        eps3=1;
        eps=[eps1, eps2,eps3];
        
        n_couches=length(eps);
        H=zeros(1,n_couches);
        
        for ih=1:n_couches-2
            H(ih+1)=h(ih);
        end
        switch pola
            case 'TE'
                p=ones(n_couches,1);
                coef_t=sqrt(eps(end)/eps(1));
            case 'TM'
                p=1./eps;
                coef_t=sqrt(eps(1)/eps(end));
        end
        
        
        
        k0=Omega/c;
        k=k0.*sqrt(eps);
        
        
        %alpha0=1/2;
        
        %alpha1=-k(1)*sin(theta);
        
        alpha1=-alpha0*sqrt(eps(1));
        
        beta0=sqrt(k.^2-alpha1.^2);
        
        for ic=2:n_couches-1
            C(:,:,ic)= [exp(-1i*beta0(ic)*H(ic)),0;0,exp(1i*beta0(ic)*H(ic))];
        end
        
        
        for it=1:n_couches-1
            s(it)=(p(it).*beta0(it)+p(it+1).*beta0(it+1))./(2*p(it).*beta0(it));
            d(it)=(p(it).*beta0(it)-p(it+1).*beta0(it+1))./(2*p(it).*beta0(it));
            T(:,:,it)=[s(it),d(it);d(it),s(it)];
        end
        
        M=T(:,:,1);
        for im=2:n_couches-1
            M=M*C(:,:,im)*T(:,:,im);
        end
        t(ilr,ili)=1./M(1,1);
        
        r(ilr,ili)=M(2,1)./M(1,1);
        
        %  if isreal (lambda0)
        %         lambg=[lambg;lambda0];
        %         fg=[fg;coef_t*abs(t(ilr,ili)).^2];
        %     end
        
    end
    
end
if norm
    out=coef_t*(abs(t)).^2;
else
    out=t;
end
%
% hold on
% figure
%
% plot(1e6*LAMB,out,'r')

