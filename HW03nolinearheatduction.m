% This is the nonlinear solver for heat equation with material nonlinearity
% clean the memory and the screen
close all
clear all; clc;

% domain
omega_l = 0.0;
omega_r = 1.0;

% -------------------------------------------------------------------------
% material properties and input data
%kappa = 1.0;

% exact solution
exact = @(x) sin(x);
exact_x = @(x) cos(x);

% nonlinear kappa
fun_kappa = @(x) 1 + x*x;
fun_dkappa = @(x) 2*x;

f = @(x) -2.0*cos(x)*cos(x)*sin(x) + sin(x)*(sin(x)*sin(x) + 1.0);
h = @(x) -1.0 * fun_kappa( exact(0) );
g = @(x) sin(1);
% -------------------------------------------------------------------------

% interpolation degree
pp = 2;

nElem_vec = [12,14,16,18,20];

cn       = length(nElem_vec);
hh_x     = zeros(cn,1);
error_l2  = zeros(cn,1);
error_h1  = zeros(cn,1);
slop_l2     = zeros(cn,1);
slop_h1     = zeros(cn,1);

for num = 1 : cn
    % number of elements
    nElem = nElem_vec(num);

    % quadrature rule
    %nqp = pp + 1;
    nqp = 6;
    [qp, wq] = Gauss(nqp, -1, 1);

    n_np = nElem * pp + 1; % number of nodal points
    n_en = pp + 1;         % number of element nodes

    IEN = zeros(n_en, nElem);

    for ee = 1 : nElem
        for aa = 1 : n_en
            IEN(aa, ee) = (ee - 1) * pp + aa;
        end
    end

    % mesh is assumbed to have uniform size hh
    hh = (omega_r - omega_l) / nElem;

    x_coor = omega_l : (hh/pp) : omega_r;

    % setup ID array based on the boundary condition
    ID = 1 : n_np;
    ID(end) = 0;

    % Setup the stiffness matrix and load vector
    % number of equations equals the number of nodes minus the number of
    % Dirichlet nodes
    n_eq = n_np - 1;

    %initial guess
    uh = [zeros(n_eq,1);g(omega_r)];

    counter = 0;
    nmax    = 20;
    error   = 1.0;

    while counter < nmax && error > 1.0e-20

        % Allocate an empty stiffness matrix and load vector
        K = sparse(n_eq, n_eq);
        F = zeros(n_eq, 1);

        % Assembly the siffness matrix and load vector
        for ee = 1 : nElem
            % Allocate zero element stiffness matrix and element load vector
            k_ele = zeros(n_en, n_en);
            f_ele = zeros(n_en, 1);

            x_ele = zeros(n_en, 1);
            d_ele = zeros(n_en, 1);
            for aa = 1 : n_en
                x_ele(aa) = x_coor( IEN(aa,ee) );
                u_ele(aa) = uh( IEN(aa,ee) );
            end

            for qua = 1 : nqp
                % geometrical mapping
                x_qua  = 0.0;
                dx_dxi = 0.0;
                u_qua  = 0.0;
                u_xi   = 0.0;

                for aa = 1 : n_en
                    x_qua  = x_qua  + x_ele(aa) * PolyBasis(pp, aa, 0, qp(qua));
                    dx_dxi = dx_dxi + x_ele(aa) * PolyBasis(pp, aa, 1, qp(qua));
                    u_qua  = u_qua  + u_ele(aa) * PolyBasis(pp, aa, 0, qp(qua));
                    u_xi   = u_xi   + u_ele(aa) * PolyBasis(pp, aa, 1, qp(qua));
                end

                dxi_dx = 1.0 / dx_dxi;

                kappa  = fun_kappa( u_qua );
                dkappa = fun_dkappa( u_qua );

                for aa = 1 : n_en
                    Na    = PolyBasis(pp, aa, 0, qp(qua));
                    Na_xi = PolyBasis(pp, aa, 1, qp(qua));
                    f_ele(aa) = f_ele(aa) + wq(qua) * Na * f(x_qua) * dx_dxi;
                    f_ele(aa) = f_ele(aa) - wq(qua) * Na_xi * kappa * u_xi * dxi_dx;
                    for bb = 1 : n_en
                        Nb    = PolyBasis(pp, bb, 0, qp(qua));
                        Nb_xi = PolyBasis(pp, bb, 1, qp(qua));
                        k_ele(aa,bb) = k_ele(aa,bb) + wq(qua) * Na_xi * kappa * Nb_xi * dxi_dx;
                        k_ele(aa,bb) = k_ele(aa,bb) + wq(qua) * Na_xi * dkappa * Nb * u_xi* dxi_dx;
                    end
                end
            end
            % end of the quadrature loop

            % distribute the entries to the global stiffness matrix and global load vector
            for aa = 1 : n_en
                LM_a = ID( IEN(aa, ee) );
                if LM_a > 0
                    F(LM_a) = F(LM_a) + f_ele(aa);
                    for bb = 1 : n_en
                        LM_b = ID( IEN(bb, ee) );
                        if LM_b > 0
                            K(LM_a, LM_b) = K(LM_a, LM_b) + k_ele(aa, bb);
                        else
                            %   x_qua =  x_coor( IEN(bb,ee) ); % obtain the Dirichlet node's physical coordinates
                            %  g_qua = g( x_qua ); % Obtain the boundary data at this point
                            %  F( LM_a ) = F( LM_a ) - k_ele(aa, bb) * g_qua;
                        end
                    end
                end
            end

            % Modify the load vector by the Natural BC
            % Note: for multi-dimensional cases, one needs to perform line or
            % surface integration for the natural BC.
            if ee == 1
                F( ID(IEN(1, ee)) ) = F( ID(IEN(1, ee)) ) + h( x_coor(IEN(1,ee)));
            end
        end

        % Solve the stiffness matrix problem
        %uh = K \ F;
        increment = K \ F;
        %[L,U,P] = lu(K);
        %y = L\(P*F);
        %uh = U \ y;
        % Append the displacement vector by the Dirichlet data
        %uh = [ uh; g(omega_r) ];

        uh = [ uh(1:end-1) + increment; g(omega_r) ];

        error = norm(F);
        counter = counter + 1;

    end
    % -------------------------------------------------------------------------
    % Now we do the postprocessing
    nqp = 6;
    [qp, wq] = Gauss(nqp, -1, 1);

    top_l2 = 0.0; bot_l2 = 0.0; top_h1 = 0.0; bot_h1 = 0.0;
    for ee = 1 : nElem
        for qua = 1 : nqp
            x_ele = zeros(n_en, 1);
            u_ele = zeros(n_en, 1);
            for aa = 1 : n_en
                x_ele(aa) = x_coor(IEN(aa, ee));
                u_ele(aa) = uh(IEN(aa, ee));
            end

            x = 0.0; dx_dxi = 0.0; duh_dxi = 0.0; uh_0 = 0.0;
            for aa = 1 : n_en
                x = x + x_ele(aa) * PolyBasis(pp, aa, 0, qp(qua));
                dx_dxi = dx_dxi + x_ele(aa) * PolyBasis(pp, aa, 1, qp(qua));
                uh_0   = uh_0 + u_ele(aa) * PolyBasis(pp, aa, 0, qp(qua));
                duh_dxi = duh_dxi + u_ele(aa) * PolyBasis(pp, aa, 1, qp(qua));
            end

            dxi_dx = 1.0 / dx_dxi;

            top_h1 = top_h1 + wq(qua) * (duh_dxi * dxi_dx - exact_x(x))^2 * dx_dxi;
            bot_h1 = bot_h1 + wq(qua) * exact_x(x)^2 * dx_dxi;

            top_l2 = top_l2 + wq(qua) * (uh_0 - exact(x))^2 * dx_dxi;
            bot_l2 = bot_l2 + wq(qua) * exact(x)^2 * dx_dxi;
        end
    end

    top_l2 = sqrt(top_l2);
    bot_l2 = sqrt(bot_l2);

    top_h1 = sqrt(top_h1);
    bot_h1 = sqrt(bot_h1);

    errorl2 = top_l2 / bot_l2;
    errorh1 = top_h1 / bot_h1;

    error_l2(num,1)  = log(errorl2);
    error_h1(num,1)  = log(errorh1);

    hh_x(num,1)    = log(hh);

    dx_hh    = zeros( 2,1);
    dy_error_l2 = zeros( 2,1);
    dy_error_h1 = zeros( 2,1);

    if (num >= 2)
        dx_hh = [hh_x(num);hh_x(num-1)];

        dy_error_l2 = [error_l2(num);error_l2(num-1)];
        dy_error_h1 = [error_h1(num);error_h1(num-1)];

        slop_l2(num) = (dy_error_l2(2) - dy_error_l2(1)) / (dx_hh(2)-dx_hh(1));
        slop_h1(num) = (dy_error_h1(2) - dy_error_h1(1)) / (dx_hh(2)-dx_hh(1));

    else
        slop_l2(num) = 0;
        slop_h1(num) = 0;
    end

end
% -------------------------------------------------------------------------
% Postprocessing plot graph
figure
%yaxis
error_x_l2 = plot(hh_x,error_l2,'--rO','LineWidth',2);
hold on
error_x_h1 = plot(hh_x,error_h1,'--rO','LineWidth',2);

%legend(error_x_l2,error_x_h1,'L_2 error convergence analysis','H_1 error convergence analysis');
xlabel('log ||hh|| ');
%ylabel('log ||e||_{l2} ');
exportgraphics(gca,['error_x_l2_pp3' '.jpg']);
T4 = table(hh_x,error_l2,error_h1,slop_l2,slop_h1,...
    'variableNames',{'hh_mesh','error_l2','error_h1','l2 convergence rate','H1 convergence rate',});
writetable(T4);
disp(T4);

% EOF

