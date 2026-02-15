% eee431 hw1 - Emir Bayer

clear; clc; close all;
rng(88);


%% part 1

n_small =100;
n_big =100000;

u_small= rand(n_small,1);
u_big =rand(n_big,1);

x_small = makexfromU(u_small);
x_big = makexfromU( u_big );

tri_pdf =@(x) (1-abs(x)).*(abs(x)<=1);

figure;
subplot(1,2,1);
histogram(x_small,'Normalization','pdf');
hold on;
fplot(tri_pdf,[-1 1],'r','LineWidth',1.3);
title('part1(b) n=100');

subplot(1,2,2);
histogram(x_big,'Normalization','pdf');
hold on;
fplot(tri_pdf,[-1 1],'r','LineWidth',1.3);
title('part1(b) n=100000');


% part 1(c)
levList =[8 16 64];

D_anl =1./(3*(levList.^2));
SQNR_anl =(levList.^2)/2;
SQNR_anl_db =10*log10(SQNR_anl);

fprintf('part 1(c)\n');
for k=1:numel(levList)
    fprintf('N=%d mse=%.4e sqnr=%.2f dB\n',levList(k),D_anl(k),SQNR_anl_db(k));
end
fprintf('\n');


% part 1(d)-(h)
xx = x_big(:);
Px_emp = mean(xx.^2);

res_unif = zeros(numel(levList),3);
res_lm = zeros(numel(levList),3);

fprintf('part 1(d)-(h)\n');
fprintf('empirical Px ~ %.4f\n',Px_emp);

for t=1:numel(levList)

    Nq = levList(t);

    dq =2/Nq;
    b_unif = linspace(-1+dq,1-dq,Nq-1);
    c_unif = linspace(-1+dq/2,1-dq/2,Nq);

    e_unif = [-inf b_unif inf];
    [~,~,iu] = histcounts(xx,e_unif);

    xq_u = c_unif(iu);
    xq_u = xq_u(:);
    Du = mean((xx-xq_u).^2);
    Sq = Px_emp/Du;
    Sq_db =10*log10(Sq);

    res_unif(t,:) =[Du Sq Sq_db];

    [b_lm,c_lm] = lloyd1d(xx,c_unif(:));

    e_lm = [-inf; b_lm(:); inf];
    [~,~,il] = histcounts(xx,e_lm');

    xq_l = c_lm(il);
    xq_l = xq_l(:);
    Dl = mean((xx-xq_l).^2);
    Sl = Px_emp/Dl;
    Sl_db =10*log10(Sl);

    res_lm(t,:) =[Dl Sl Sl_db];

    % pmf plots only for N=16 (representative for report)
    if Nq==16
        pmf_u= histcounts(iu,1:Nq+1,'Normalization','probability');
        figure;
        stem(c_unif,pmf_u,'filled');
        title(sprintf('pmf uniform N=%d',Nq));

        pmf_l= histcounts(il,1:Nq+1,'Normalization','probability');
        figure;
        stem(c_lm,pmf_l,'filled');
        title(sprintf('pmf lloyd-max N=%d',Nq));
    end

end

fprintf('\nuniform quant:\n');
for t=1:numel(levList)
    fprintf('N=%d mse=%.6f sqnr=%.2f dB\n',levList(t),res_unif(t,1),res_unif(t,3));
end

fprintf('\nlloyd-max quant:\n');
for t=1:numel(levList)
    fprintf('N=%d mse=%.6f sqnr=%.2f dB\n',levList(t),res_lm(t,1),res_lm(t,3));
end


%% part 2

picList ={"flower.jpg","foliage.tif"};
Ngray =[2 4 8 16 32];
polyOrd =4;

for p=1:numel(picList)

    I0 = imread(picList{p});
    if ndims(I0)==3
        I0 = rgb2gray(I0);
    end
    I0 = double(I0);

    [h,w] = size(I0);
    Pxv = var(I0(:));

    g=0:255;
    cnt= histcounts(I0(:),[g 256]);
    pmf = cnt/sum(cnt);

    pc= polyfit(g,pmf,polyOrd);
    pdf_hat= polyval(pc,g);
    pdf_hat(pdf_hat<1e-6)=1e-6;
    pdf_hat = pdf_hat/sum(pdf_hat);

    figure;
    subplot(2,1,1);
    bar(g,pmf);
    title(sprintf('%s histogram',picList{p}));

    subplot(2,1,2);
    plot(g,pdf_hat,'LineWidth',1.3);
    title(sprintf('%s poly pdf',picList{p}));

    Dm=zeros(size(Ngray));
    Sm=zeros(size(Ngray));
    Iq= cell(size(Ngray));

    fprintf('\nimage: %s\n',picList{p});

    figure;
    for j=1:numel(Ngray)

        Ng = Ngray(j);

        step=256/Ng;
        qind =floor(I0/step);
        qind(qind>=Ng)=Ng-1;

        IQ = qind*step + step/2;
        Iq{j}=IQ;

        Dm(j)= mean((I0(:)-IQ(:)).^2);
        Sm(j)= 10*log10(Pxv/Dm(j));

        subplot(2,3,j);
        imshow(uint8(IQ));
        title(sprintf('%s N=%d',picList{p},Ng));
    end

    figure;
    plot(Ngray,Sm,'-o','LineWidth',1.2);
    xlabel('N'); ylabel('sqnr');
    title(sprintf('%s sqnr vs N',picList{p}));

    % fft
    F0 = fftshift(fft2(I0));
    S0 = log(1+abs(F0));

    figure;
    subplot(2,3,1);
    imagesc(S0); colormap gray; axis image off;
    title('fft orig');

    for j=1:numel(Ngray)
        Fq = fftshift(fft2(Iq{j}));
        Sq = log(1+abs(Fq));
        subplot(2,3,j+1);
        imagesc(Sq); colormap gray; axis image off;
        title(sprintf('fft N=%d',Ngray(j)));
    end

    % LPF for N=8
    idx8 = find(Ngray==8,1);
    I8 = Iq{idx8};

    [U,V] = meshgrid(-w/2:w/2-1, -h/2:h/2-1);
    Hlp = (abs(U)<0.1*w) & (abs(V)<0.1*h);

    F0f = fftshift(fft2(I0));
    F8f = fftshift(fft2(I8));

    F0f = F0f.*Hlp;
    F8f = F8f.*Hlp;

    I0lp = real(ifft2(ifftshift(F0f)));
    I8lp = real(ifft2(ifftshift(F8f)));

    figure;
    subplot(1,3,1); imshow(uint8(I0)); title('orig');
    subplot(1,3,2); imshow(uint8(I0lp)); title('lpf orig');
    subplot(1,3,3); imshow(uint8(I8lp)); title('lpf N=8');

end


%% helpers

function x = makexfromU(u)
    x = zeros(size(u));
    L = (u<=0.5);
    R = ~L;
    x(L)= sqrt(2*u(L)) -1;
    x(R)= 1 - sqrt(2*(1-u(R)));
end

function [b,c] = lloyd1d(x,init)
    N = length(init);
    c = sort(init(:));
    tol=1e-6;

    for it=1:100
        b = (c(1:N-1)+c(2:N))/2;
        e=[-inf; b; inf];
        [~,~,idx]= histcounts(x,e');

        nc=c;
        for k=1:N
            vals = x(idx==k);
            if ~isempty(vals)
                nc(k)= mean(vals);
            end
        end
        if norm(nc-c)/max(norm(c),eps) < tol
            c=nc; break;
        end
        c=nc;
    end
    b= (c(1:N-1)+c(2:N))/2;
end
