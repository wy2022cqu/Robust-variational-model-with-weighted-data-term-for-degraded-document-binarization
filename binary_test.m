function u = binary_test(u0,c_u,r,l,lmada,c_W,alpha)

    % fixed parameter 
    iterNum = 30;
    delta = 0.25;
    bp = -0.5;%(bp=-a,a=0.5 in our paper) 

    % obtain D, W0, sigma
    lc = local_contrast(u0,l);
    [uf,ub,c] = generate_hard(u0,r);
    mad = MAD(u0);
    W0 = (abs(ub-uf)).*lc;
    

    % initialization
    u = u0;
    W = lc;


    for t1 = 1:iterNum

        % regularization
        old_u = u;
        u_regular = PM(u,mad,1);

        % F(u) and -F_u(u)
        [F_u,B_u] = F(u, bp, c);

        % evolution
        W = W + delta* (-c_W*F_u + lmada*(W0 - W));
        W = max(0,W);
        u = u + delta* (alpha*u_regular + c_u* W.* B_u);

        % rela_error
        rela_error(t1) = norm(u-old_u,'fro')/norm(old_u,'fro');
        if rela_error(t1) <= 0.0001
            t1
            break 
        end
    end

end

%% 函数区域

function lc = local_contrast(u0,l)

    % ---- compensated image ----
    B = convMax(single(u0),l);
    B = guidedfilter(u0, B, 20, 0.0001);
    C = u0 - B;
    C = C - min(C(:));

    % % ---- local contrast ----
    Max = convMax(single(C),l);
    Min = convMax(single(1-C),l);
    Min = 1-Min;
    lc = (Max - Min)./(Max + Min + 10e-3);
end

function [uf,ub,c] = generate_hard(u,r)
    G = fspecial('gaussian', round(2*r)*2+1, r);
    us = imfilter(u, G, 'replicate');

    maskf = (u < us - 1e-5);
    maskb = ~maskf;

    maskf_d = imfilter(single(maskf), G, 'replicate');
    maskb_d = imfilter(single(maskb), G, 'replicate');
    
    % local text/background mean
    uf = imfilter(u .* single(maskf), G, 'replicate') ./ (maskf_d + 1e-12);
    ub = imfilter(u .* single(maskb), G, 'replicate') ./ (maskb_d + 1e-12);

    %local threshold
    c = 0.5 * (ub + uf);
end



function out = PM(u,mad,c)

    sigma = c*mad;

    u_padded = padarray(u, [1 1], 'symmetric');

    [dx,dy] = gradient(u_padded);
    grad = hypot(dx, dy);   

    g = ones(size(grad), 'like', grad);
    mask = grad > sigma;
    g(mask) = sigma ./ (grad(mask) + 1e-6);

    [dxx,~] = gradient(g .* dx);
    [~,dyy] = gradient(g .* dy);

    out = dxx(2:end-1, 2:end-1) + dyy(2:end-1, 2:end-1);

end


function [F_u,B_u] = F(u, a, T)
    % a = bp = -0.5  

    c1 = (T - a) / pi + a * log(abs(a)) - a;
    c2 = (T - a) / pi + (2*T - a) + a * log(abs(a));

    mask1 = (u <  a);
    mask2 = (u >= a) & (u <= 2*T - a);
    mask3 = (u >  2*T - a);

    epsl = 1e-12;
    u_safe         = max(abs(u), epsl);           % 用于 ln(u)
    twoT_minus_u   = 2*T - u;
    twoT_minus_u_s = max(abs(twoT_minus_u), epsl);% 用于 ln(2T-u)

    F1 = -a.*log(u_safe) + u + c1;
    F2 = -((T - a)/pi) .* cos((pi.*(u - T))./ (T - a));
    F3 = -a.*log(twoT_minus_u_s) - u + c2;

    % F(u) (a = bp = -0.5  )
    F_u = -(F1.*mask1 + F2.*mask2 + F3.*mask3);

    % -F_u(u) (a = bp = -0.5 )
    up = 2*T - a; omega = pi./(T - a);
    u_b = a./(-(u.*double(u<=a))- a.*double(u>a)) + 1;
    u_01 = sin(omega.*(u-T)).*double(u<=up).*double(u>=a);
    u_u = a./(2*T - u.*double(u>=up) + (-2*T+a).*double(u<up))-1;
    B_u = u_b + u_01 + u_u;

 
end

function mad = MAD(u)
    [ux,uy] = gradient(u);
    grad = sqrt(ux.^2+uy.^2);

    temp = abs(grad - median(grad(:)));
    mad = 1.4826*median(temp(:));

end