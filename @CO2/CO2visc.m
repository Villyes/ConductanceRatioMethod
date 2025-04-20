function eta=CO2visc(rho,T)
    persistent a epsByK
    if isempty(a)
        vars=load('ViscosityConstants.mat','a','e');
        a=vars.a;
        e=vars.e;
        
        epsByK=251.196;
    end
    
    sz=implExpSize(rho,T);
    [rho,T]=normalize(sz,rho,T);
    
    ecs=exp(sum(a.*log(T./epsByK).^(0:4)',1));
    eta_0=1.00697.*sqrt(T)./ecs;
    
    DeltaEps_g=sum(e.*rho.^(1:7)',1);
    
end



function sz=implExpSize(varargin)
    dims=max(cellfun(@(x) ndims(x),varargin));
    szCell=cellfun(@(x) padsz(dims,x),varargin,'UniformOutput',false);
    sz=cell2mat(szCell');

    %all dimensions not equal to 1 must be the same
    not1=sz~=1;
    check=arrayfun(@(x) all(sz(not1(:,x),x)==sz(find(not1(:,x),1),x)),1:size(sz,2));

    if any(~check)
        throw(MException('incompatibleArrays','Arrays have incompatible sizes'));
    else
        sz=max(sz,[],1);
    end
end
        
        
function varargout=normalize(sz,varargin)
    n=prod(sz);
    pad=@(x) padsz(numel(sz),x);

    varargout=cellfun(@(x) reshape(repmat(x,sz./pad(x)),1,n),varargin,'UniformOutput',false);
end


function padsz=padsz(dims,x)
    padsz=[size(x),ones(1,dims-ndims(x))];
end




