clearvars
global lx ly lz bta rho c fd S V

lx = 7; ly = 5; lz = 3;
rho = 1.2;
c = 340;
bta = (rho*c)/(6e2*1820);
S = 2*lx*ly + 2*lx*lz + 2*ly*lz;
V = lx*ly*lz;
fd = oct_fraction(20, 150, 12);

lix = 2; uix = 7;
liy = 2; uiy = 5;
liz = 2; uiz = 3;

N = 10;
limit = 10;
a = 0.5;
max_iter = 60;

dane = zeros(N, 6);
Lp = zeros(N, 1);
trial = zeros(N, 1);

% === HISTORIA LICZNIKA NIEPOWODZEŃ ===
trial_history = zeros(N, max_iter);   

for i = 1:N
    dane(i,1:6) = [ ...
        lix + rand*(uix - lix), ...
        liy + rand*(uiy - liy), ...
        liz + rand*(uiz - liz), ...
        lix + rand*(uix - lix), ...
        liy + rand*(uiy - liy), ...
        liz + rand*(uiz - liz)];
    [stdLp,~,~] = room_student(dane(i,:));
    Lp(i) = stdLp;
end

% === GLOBALNE NAJLEPSZE ROZWIĄZANIE ===
[global_best_val, idx] = min(Lp);     %%% FIXED
global_best_dane = dane(idx,:);       %%% FIXED

for iter = 1:max_iter
    % FAZA 1
    for i = 1:N
        k = randi(N);
        while k == i
            k = randi(N);
        end
        phi = (2*rand(1,6)-1)*a;
        vi = dane(i,:) + phi .* (dane(i,:) - dane(k,:));
        vi = max(min(vi, [uix uiy uiz uix uiy uiz]), [lix liy liz lix liy liz]);
        [new_stdLp,~,~] = room_student(vi);

        if new_stdLp < Lp(i)
            dane(i,:) = vi;
            Lp(i) = new_stdLp;
            trial(i) = 0;
        else
            trial(i) = trial(i) + 1;
        end
    end

    % FAZA 2
    fit = exp(-Lp / mean(Lp));
    p = fit / sum(fit);
    for i = 1:N
        r = rand; acc = 0;
        for k = 1:N
            acc = acc + p(k);
            if r <= acc, 
                break; 
            end
        end
        phi = (2*rand(1,6)-1)*a;
        vi = dane(i,:) + phi .* (dane(i,:) - dane(k,:));
        vi = max(min(vi, [uix uiy uiz uix uiy uiz]), [lix liy liz lix liy liz]);
        [new_stdLp,~,~] = room_student(vi);
        if new_stdLp < Lp(i)
            dane(i,:) = vi;
            Lp(i) = new_stdLp;
            trial(i) = 0;
        else
            trial(i) = trial(i) + 1;
        end
    end

    % FAZA 3
    for i = 1:N
        if trial(i) > limit
            dane(i,1:6) = [ ...
                lix + rand*(uix - lix), ...
                liy + rand*(uiy - liy), ...
                liz + rand*(uiz - liz), ...
                lix + rand*(uix - lix), ...
                liy + rand*(uiy - liy), ...
                liz + rand*(uiz - liz)];
            [stdLp,~,~] = room_student(dane(i,:));
            Lp(i) = stdLp;
            trial(i) = 0;
        end
    end

    % === AKTUALIZACJA HISTORII LICZNIKÓW ===
    trial_history(:, iter) = trial;   %%% FIXED

    % === ŚLEDZENIE NAJLEPSZEGO ROZWIĄZANIA ===
    [iter_best_val, iter_best_idx] = min(Lp);
    if iter_best_val < global_best_val
        global_best_val = iter_best_val;
        global_best_dane = dane(iter_best_idx,:);
    end

    fprintf("Iteracja %d: Najlepszy stdLp = %.4f\n", iter, global_best_val);  %%% FIXED
end

fprintf("\nNajlepsze rozwiązanie (stdLp = %.4f):\n", global_best_val);
disp(global_best_dane);

% === WYKRES CHARAKTERYSTYKI AMPLITUDOWO-CZĘSTOTLIWOŚCIOWEJ ===
r0 = global_best_dane(1:3);
r = global_best_dane(4:6);
GF = GreenFunction_OK(r0, r);
Lp = 20 * log10(abs(1.21 * GF) / 2e-5);

figure;
semilogx(fd, Lp, '-o', 'LineWidth', 1.5);
xlabel('Częstotliwość [Hz]');
ylabel('Poziom ciśnienia akustycznego L_p [dB]');
title('Charakterystyka amplitudowo-częstotliwościowa najlepszego układu');
grid on;

% === WYKRES LICZNIKÓW NIEPOWODZEŃ ===
figure;
hold on;
for i = 1:N
    plot(1:max_iter, trial_history(i,:), 'DisplayName', sprintf('Pszczoła %d', i));
end
hold off;
xlabel('Numer iteracji');
ylabel('Licznik niepowodzeń');
legend show;
title('Ewolucja liczników niepowodzeń pszczół');
grid on;
