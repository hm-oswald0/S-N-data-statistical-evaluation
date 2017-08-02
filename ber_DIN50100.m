function [N,d_sigma,k,s_log_N] = ber_DIN50100(N_A,L_a,n_stdv,k_s)
%% Prüfung nach Horizonten- oder Perlenschnurverfahren
    if iscell(N_A)
        N_A = cell2mat(N_A);
    end
    if iscell(L_a)
        L_a = cell2mat(L_a);
    end
    test_var = unique(L_a);
    if length(test_var) == 2
        var = 1;
    else
        var = 2;
    end
n = length(N_A);
%% Berechnung nach Perlenschnurverfahren
if var == 2
    x   = log10(L_a);
    y   = log10(N_A);
    xy = x.*y;
    xx = x.*x;
    % Berechnungssummen
    sx = sum(x);
    sy = sum(y);
    sxy = sum(xy);
    sxx = sum(xx);
    % Berechnungskonstanten
    a1  = (n*sxy-sx*sy) / (n*sxx-sx^2);
    % Neigung und Lage
    if exist('k_s')
        if  isempty(k_s) == 0
            k = k_s;
        else
            k = -a1;
        end
    else
        k = -a1;
    end
    if exist ('n_stdv')
        if  isempty(n_stdv) == 0
            n_stdv = n_stdv;
        else
            n_stdv = 2;
        end
    else
        n_stdv = 2;
    end
    a0  = 1/n * (sy+k*sx);
    C = 10^a0;
    % 50%-Linie
    N_50(1) = C*1.^-k;
    N_50(2) = C*1000.^-k;
    L_a_fiktiv = 1;
    N_i_fiktiv = N_A.*(L_a_fiktiv./L_a).^-k;
    log_N_i_fiktiv = log10(N_i_fiktiv);
    log_N_50_fiktiv = 1/n * sum(log_N_i_fiktiv);
    N_50_fiktiv = 10^log_N_50_fiktiv;
    % Standardabweichung
    s_log_N = sqrt(1/(n-2) * sum((log_N_i_fiktiv-log_N_50_fiktiv).^2));
    s_log_N_korr = s_log_N * (n-1.74)/(n-2);
    % Streuspanne
    T_N = 10^(4*s_log_N_korr);
    % 10%-Linie
    N_10_fiktiv = 10^(log_N_50_fiktiv-n_stdv*s_log_N_korr);
    N_10(1) = N_10_fiktiv * (1/L_a_fiktiv).^-k;
    N_10(2) = N_10_fiktiv * (1000/L_a_fiktiv).^-k;
    % 90%-Linie
    N_90_fiktiv = 10^(log_N_50_fiktiv+n_stdv*s_log_N_korr);
    N_90(1) = N_90_fiktiv * (1/L_a_fiktiv).^-k;
    N_90(2) = N_90_fiktiv * (1000/L_a_fiktiv).^-k;
    % Übergabe der Daten
    k_end = [k;k;k];
    s_log_N_korr_end = [s_log_N_korr;s_log_N_korr];
    T_N_end = [T_N;T_N];
    N_end_pre = [N_10;N_50;N_90];
    L_a_end_pre = [1,1000;1,1000;1,1000];
%% Berechnung nach Horizontenverfahren
elseif var == 1
    if exist ('n_stdv')
        if  isempty(n_stdv) == 0
            n_stdv = n_stdv;
        else
            n_stdv = 2;
        end
    else
        n_stdv = 2;
    end
    L_a1 = L_a(L_a == max(L_a));
    L_a2 = L_a(L_a == min(L_a));
    n1 = length(L_a1);
    n2 = length(L_a2);
    N_A1 = N_A(L_a == max(L_a));
    N_A2 = N_A(L_a == min(L_a));
    % Lasthorizont 1
    log_N_50_L_a1 = 1/n1 * sum(log10(N_A1));
    N_50_L_a1 = 10^log_N_50_L_a1;
    s_log_N_L_a1 = sqrt(1/(n1-1) * sum((log10(N_A1)-log_N_50_L_a1).^2));
    s_log_N_L_a1_korr = s_log_N_L_a1 * (n1-0.74)/(n1-1);
    T_N_L_a1 = 10^(4*s_log_N_L_a1_korr);
    % Lasthorizont 2
    log_N_50_L_a2 = 1/n2 * sum(log10(N_A2));
    N_50_L_a2 = 10^log_N_50_L_a2;
    s_log_N_L_a2 = sqrt(1/(n2-1) * sum((log10(N_A2)-log_N_50_L_a2).^2));
    s_log_N_L_a2_korr = s_log_N_L_a2 * (n2-0.74)/(n2-1);
    T_N_L_a2 = 10^(4*s_log_N_L_a2_korr);
    % Werte 50%
    N_50(1) = N_50_L_a1;
    N_50(2) = N_50_L_a2;
    La1 = unique(L_a1);
    La2 = unique(L_a2);
    k_50 = -log10(N_50(1)/N_50(2))/log10(La1/La2);
    % Werte 10%
    N_10(1) = 10^(log_N_50_L_a1-n_stdv*s_log_N_L_a1_korr);
    N_10(2) = 10^(log_N_50_L_a2-n_stdv*s_log_N_L_a2_korr);
    k_10 = -log10(N_10(1)/N_10(2))/log10(La1/La2);
    % Werte 90%
    N_90(1) = 10^(log_N_50_L_a1+n_stdv*s_log_N_L_a1_korr);
    N_90(2) = 10^(log_N_50_L_a2+n_stdv*s_log_N_L_a2_korr);
    k_90 = -log10(N_90(1)/N_90(2))/log10(La1/La2);
    % Übergabe der Daten
    k_end = [k_10;k_50;k_90];
    s_log_N_korr_end = [s_log_N_L_a1_korr;s_log_N_L_a2_korr];
    T_N_end = [T_N_L_a1;T_N_L_a2];
    N_end_pre = [N_10;N_50;N_90];
    L_a_end_pre = [La1,La2;La1,La2;La1,La2];
end
% Datenumrechnung zum Plot
for ii = 1:1:3
    N = N_end_pre(ii,1);
    L = L_a_end_pre(ii,1);
    k(ii) = k_end(ii);
    Nu = 10^4;
    No = 5*10^6;
    N_end(ii,:) = [Nu,No];
    L_a_end(ii,:) = [L*(Nu/N)^(-1/k(ii));L*(No/N)^(-1/k(ii))];
end
N = [10^4; 5*10^6];
d_sigma = [L_a_end(1,1), L_a_end(2,1), L_a_end(3,1);...
           L_a_end(1,2), L_a_end(2,2), L_a_end(3,2)];
k = k_end;
s_log_N = s_log_N_korr_end;
end