

clear,clc;

k = 4;
vm = 25;
vti = 100;

% Nseq = 100;
NI = 31;

alphas = -1:0.1:1;
betas = [0, .1, .2, -.1, -.2];

%% Figure 1

MUvarAn = [];

for exp = 1:100
    
    M = gamrnd(k,sqrt(vm/k),[NI 1]);
    T = gamrnd(k,sqrt(vti/k),[NI-1 1]);
    
    metronome = mean(T); % Molly selected mean(gamma) in general which is ~ 14. ??
    
    beta = 0;
    
    for a = 1:length(alphas)
        
        alpha = alphas(a);
        
        A = [];
        I = [];
        % first tap -----------------------------------------
        n = 1;
        I(n) = T(n) + M(n+1) - M(n); % no asynchrony at first ITI ? ((Molly selected I(1) = M(1)))
        A(n) = I(n) - metronome;
        
        % second-to-last taps ------------------------------
        for n = 2:30
            A(n) = I(n-1) - metronome;
            I(n) = T(n) - alpha*A(n) - beta*A(n-1) + M(n+1) - M(n);     % tap to correct previous tap's error
        end
        
        I = I';
        A = A';
        
        allexp(exp).params(a).alpha = alpha;
        allexp(exp).params(a).beta = beta;
        allexp(exp).params(a).varAn = var(A);
        
    end
 
    MUvarAn = horzcat(MUvarAn,[allexp(exp).params.varAn].');
   
end

MUvarAn = mean(MUvarAn,2); % if A is a matrix, then mean(A,2) is a column vector containing the mean of each row.


% subplot(2,3,1);
EK_quickplot(alphas,MUvarAn,'alpha','var(An)','mean from 100 exp');%hold on;
set(gca, 'YLim', [0 1400]); set(gca, 'XLim', [min(alphas) max(alphas)]);


% subplot(2,3,2)
% EK_quickplot([params.alpha],[params.ACorr1],'alpha','p(An,An+1)',[]);hold on;
% set(gca, 'YLim', [-1 1]); set(gca, 'XLim', [0 2]);
% subplot(2,3,3)
% EK_quickplot([params.alpha],[params.ACorr2],'alpha','p(An,An+2)',[]);hold on;
% set(gca, 'YLim', [-1 1]); set(gca, 'XLim', [0 2]);
% subplot(2,3,4)
% EK_quickplot([params.alpha],[params.varIn],'alpha','var(In)',[]);hold on;
% set(gca, 'YLim', [0 1400]); set(gca, 'XLim', [0 2]);
% subplot(2,3,5)
% EK_quickplot([params.alpha],[params.ICorr1],'alpha','p(In,In+1)',[]);hold on;
% set(gca, 'YLim', [-1 1]); set(gca, 'XLim', [0 2]);
% subplot(2,3,6)
% EK_quickplot([params.alpha],[params.ICorr2],'alpha','p(In,In+2)',[]);hold on;
% set(gca, 'YLim', [-1 1]); set(gca, 'XLim', [0 2]);

% end





%
%     params(exp).ACorr1 = corr(A(1:end-1), A(2:end));
%     params(exp).ACorr2 = corr(A(1:end-2), A(3:end));
%
%     params(exp).varIn = var(I);
%     params(exp).ICorr1 = corr(I(1:end-1), I(2:end));
%     params(exp).ICorr2 = corr(I(1:end-2), I(3:end));
%


