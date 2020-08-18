

clear,clc;

k = 4;
vm = 25;
vti = 100;

NI = 31;

alphas = 0.2:0.1:1.8;
betas = [0, .1, .2, -.1, -.2];

params = struct();

for exp = 1:100
    
    M = gamrnd(k,sqrt(vm/k),[NI 1]);
    T = gamrnd(k,sqrt(vti/k),[NI-1 1]);
    metronome = mean(T); 

    for b = 1:length(betas)
        beta = betas(b);
        for a = 1:length(alphas)
            alpha = alphas(a);
            A = [];
            I = [];
            for n = 1:30
                if n == 1
%                     A(n) = gamrnd(k,sqrt(vm/k),[1 1]);
                    A(n) = 0;
                    Tx(n) = T(n) -  alpha * A(n);
                    I(n) = Tx(n) + M(n+1) - M(n);
                else
                    A(n) = sum(I(1:n-1)) - metronome*(n-1);
                    Tx(n) = T(n)  - alpha * A(n) - beta * A(n-1);
                    I(n) = Tx(n) + M(n+1) - M(n);
                end
            end
            I = I'; A = A';
            
            params.alpha(a) = alpha;
            
            params.varAn(exp,b,a) = var(A);
            params.ACorr1(exp,b,a) = corr(A(1:end-1), A(2:end));
            params.ACorr2(exp,b,a) = corr(A(1:end-2), A(3:end));
            
            params.varIn(exp,b,a) = var(I);
            params.ICorr1(exp,b,a) = corr(I(1:end-1), I(2:end));
            params.ICorr2(exp,b,a) = corr(I(1:end-2), I(3:end));

        end % a
    end % b
end % exp

mu_varAn = mean(params.varAn,1); mu_varAn = squeeze(mu_varAn(1, :, :));
mu_ACorr1 = mean(params.ACorr1,1); mu_ACorr1 = squeeze(mu_ACorr1(1, :, :));
mu_ACorr2 = mean(params.ACorr2,1); mu_ACorr2 = squeeze(mu_ACorr2(1, :, :));
mu_varIn = mean(params.varIn,1); mu_varIn = squeeze(mu_varIn(1, :, :));
mu_ICorr1 = mean(params.ICorr1,1); mu_ICorr1 = squeeze(mu_ICorr1(1, :, :));
mu_ICorr2 = mean(params.ICorr2,1); mu_ICorr2 = squeeze(mu_ICorr2(1, :, :));


figure
for b = 1:length(betas)
    
    subplot(2,3,1);
    EK_quickplot(params.alpha,mu_varAn(b,:),'alpha','var(An)',[]); hold on;
    set(gca, 'YLim', [0 1400]); set(gca, 'XLim', [0 2]);
    
    subplot(2,3,2)
    EK_quickplot(params.alpha, mu_ACorr1(b,:),'alpha','p(An,An+1)',[]);hold on;
    set(gca, 'YLim', [-1 1]); set(gca, 'XLim', [0 2]);  yline(0);
    
    subplot(2,3,3)
    EK_quickplot(params.alpha, mu_ACorr2(b,:),'alpha','p(An,An+2)',[]);hold on;
    set(gca, 'YLim', [-1 1]); set(gca, 'XLim', [0 2]);  yline(0);
    
    subplot(2,3,4)
    EK_quickplot(params.alpha, mu_varIn(b,:),'alpha','var(In)',[]);hold on;
    set(gca, 'YLim', [0 1400]); set(gca, 'XLim', [0 2]);  
    
    subplot(2,3,5)
    EK_quickplot(params.alpha, mu_ICorr1(b,:),'alpha','p(In,In+1)',[]);hold on;
    set(gca, 'YLim', [-1 1]); set(gca, 'XLim', [0 2]);  yline(0);
    
    subplot(2,3,6)
    EK_quickplot(params.alpha, mu_ICorr2(b,:),'alpha','p(In,In+2)',[]);hold on;
    set(gca, 'YLim', [-1 1]); set(gca, 'XLim', [0 2]);  yline(0);
    
    
end


