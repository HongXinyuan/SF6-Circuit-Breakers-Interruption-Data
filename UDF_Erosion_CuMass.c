// Define mass source term in x direction
DEFINE_SOURCE(mass_vx_source, c, t, dS, eqn)
{
    real x[ND_ND];  // Coordinate array, used to store the grid point coordinates
    C_CENTROID(x, c, t);  // Get cell centroid coordinates
    real x_lim = x[0];
    real y_lim = x[1];
    real z_lim = x[2];
    real pvap;
    real T_vap = 2835.0;  // Define evaporation temperature
    real L_vap = 3.6e5;  // Define latent heat
    real M_Cu = 0.063546;  // Define molecular weight of copper
    real Pcu = C_P(c, t);  // Get cell pressure
    real T = C_T(c, t);  // Get cell temperature
    real source = 0.0;  // Initialize source term
    real delta_t = CURRENT_TIMESTEP;  // Get current time step size
    real nx = 0.0, ny = 0.0, nz = 0.0;
    pvap = exp(L_vap * (1.0 / T_vap - 1.0 / T) / 56.94);  // Calculate vapor pressure
    real omiga = (pvap * M_Cu) / (pvap * M_Cu + (1.0 - pvap) * 0.14606);  // Calculate mass fraction parameter
    real b = RP_Get_Real("flow-time");  // Get current flow time

    if (b <= 0.000125)
    {
        source = 0.0;  // If time is less than 0.000125, source term is 0
    }
    else {
        // Check whether in the top-surface ablation region
        if (is_on_top_surface(x_lim, y_lim, z_lim))
        {
            // Get the normal direction of the top surface
            get_top_surface_normal(x_lim, y_lim, z_lim, &nx, &ny, &nz);

            if (T >= T_vap)
            {
                source = (0.06 * 4.0 * 3.1415926 * net(T) / 4.6e6) * nx;  // x-direction component
            }
            else
            {
                source = (0.06 * 4.0 * 3.1415926 * net(T) / 4.6e6 - omiga * Pcu * M_Cu / (8.314 * T * delta_t)) * nx;  // x-direction component
            }
        }
        // Check whether in the side-surface ablation region
        else if (is_on_side_surface(x_lim, y_lim, z_lim))
        {
            // Get the normal direction of the side surface
            get_side_surface_normal(x_lim, y_lim, z_lim, &nx, &ny, &nz);

            if (T >= T_vap)
            {
                source = (0.06 * 4.0 * 3.1415926 * net(T) / 4.6e6) * nx;  // x-direction component
            }
            else
            {
                source = (0.06 * 4.0 * 3.1415926 * net(T) / 4.6e6 - omiga * Pcu * M_Cu / (8.314 * T * delta_t)) * nx;  // x-direction component
            }
        }
        else
        {
            source = 0.0;  // If not in the ablation region, source term is 0
        }
    }

    dS[eqn] = 0.0;  // Set derivative of source term to 0
    return source;  // Return source term value
}

// Define mass source term in y direction
DEFINE_SOURCE(mass_vy_source, c, t, dS, eqn)
{
    real x[ND_ND];
    C_CENTROID(x, c, t);  // Get cell centroid coordinates
    real x_lim = x[0];
    real y_lim = x[1];
    real z_lim = x[2];
    real pvap;
    real T_vap = 2835.0;  // Define evaporation temperature
    real L_vap = 3.6e5;  // Define latent heat
    real M_Cu = 0.063546;  // Define molecular weight of copper
    real Pcu = C_P(c, t);  // Get cell pressure
    real T = C_T(c, t);  // Get cell temperature
    real source = 0.0;  // Initialize source term
    real delta_t = CURRENT_TIMESTEP;  // Get current time step size
    real nx = 0.0, ny = 0.0, nz = 0.0;
    pvap = exp(L_vap * (1.0 / T_vap - 1.0 / T) / 56.94);  // Calculate vapor pressure
    real omiga = (pvap * M_Cu) / (pvap * M_Cu + (1.0 - pvap) * 0.14606);  // Calculate mass fraction parameter
    real b = RP_Get_Real("flow-time");  // Get current flow time

    if (b <= 0.000125)
    {
        source = 0.0;  // If time is less than 0.000125, source term is 0
    }
    else {
        // Check whether in the top-surface ablation region
        if (is_on_top_surface(x_lim, y_lim, z_lim))
        {
            // Get the normal direction of the top surface
            get_top_surface_normal(x_lim, y_lim, z_lim, &nx, &ny, &nz);

            if (T >= T_vap)
            {
                source = (0.06 * 4.0 * 3.1415926 * net(T) / 4.6e6) * ny;  // y-direction component
            }
            else
            {
                source = (0.06 * 4.0 * 3.1415926 * net(T) / 4.6e6 - omiga * Pcu * M_Cu / (8.314 * T * delta_t)) * ny;  // y-direction component
            }
        }
        // Check whether in the side-surface ablation region
        else if (is_on_side_surface(x_lim, y_lim, z_lim))
        {
            // Get the normal direction of the side surface
            get_side_surface_normal(x_lim, y_lim, z_lim, &nx, &ny, &nz);

            if (T >= T_vap)
            {
                source = (0.06 * 4.0 * 3.1415926 * net(T) / 4.6e6) * ny;  // y-direction component
            }
            else
            {
                source = (0.06 * 4.0 * 3.1415926 * net(T) / 4.6e6 - omiga * Pcu * M_Cu / (8.314 * T * delta_t)) * ny;  // y-direction component
            }
        }
        else
        {
            source = 0.0;  // If not in the ablation region, source term is 0
        }
    }

    dS[eqn] = 0.0;  // Set derivative of source term to 0
    return source;  // Return source term value
}

// Define mass source term in z direction
DEFINE_SOURCE(mass_vz_source, c, t, dS, eqn)
{
    real x[ND_ND];
    C_CENTROID(x, c, t);  // Get cell centroid coordinates
    real x_lim = x[0];
    real y_lim = x[1];
    real z_lim = x[2];
    real pvap;
    real T_vap = 2835.0;  // Define evaporation temperature
    real L_vap = 3.6e5;  // Define latent heat
    real M_Cu = 0.063546;  // Define molecular weight of copper
    real Pcu = C_P(c, t);  // Get cell pressure
    real T = C_T(c, t);  // Get cell temperature
    real source = 0.0;  // Initialize source term
    real delta_t = CURRENT_TIMESTEP;  // Get current time step size
    real nx = 0.0, ny = 0.0, nz = 0.0;
    pvap = exp(L_vap * (1.0 / T_vap - 1.0 / T) / 56.94);  // Calculate vapor pressure
    real omiga = (pvap * M_Cu) / (pvap * M_Cu + (1.0 - pvap) * 0.14606);  // Calculate mass fraction parameter
    real b = RP_Get_Real("flow-time");  // Get current flow time

    if (b <= 0.000125)
    {
        source = 0.0;  // If time is less than 0.000125, source term is 0
    }
    else {
        // Check whether in the top-surface ablation region
        if (is_on_top_surface(x_lim, y_lim, z_lim))
        {
            // Get the normal direction of the top surface
            get_top_surface_normal(x_lim, y_lim, z_lim, &nx, &ny, &nz);

            if (T >= T_vap)
            {
                source = (0.06 * 4.0 * 3.1415926 * net(T) / 4.6e6) * nz;  // z-direction component
            }
            else
            {
                source = (0.06 * 4.0 * 3.1415926 * net(T) / 4.6e6 - omiga * Pcu * M_Cu / (8.314 * T * delta_t)) * nz;  // z-direction component
            }
        }
        // Check whether in the side-surface ablation region
        else if (is_on_side_surface(x_lim, y_lim, z_lim))
        {
            // Get the normal direction of the side surface
            get_side_surface_normal(x_lim, y_lim, z_lim, &nx, &ny, &nz);

            if (T >= T_vap)
            {
                source = (0.06 * 4.0 * 3.1415926 * net(T) / 4.6e6) * nz;  // z-direction component
            }
            else
            {
                source = (0.06 * 4.0 * 3.1415926 * net(T) / 4.6e6 - omiga * Pcu * M_Cu / (8.314 * T * delta_t)) * nz;  // z-direction component
            }
        }
        else
        {
            source = 0.0;  // If not in the ablation region, source term is 0
        }
    }

    dS[eqn] = 0.0;  // Set derivative of source term to 0
    return source;  // Return source term value
}

// Define mass source term
DEFINE_SOURCE(mass1, c, t, dS, eqn)
{
    real x[ND_ND];
    C_CENTROID(x, c, t);  // Get cell centroid coordinates
    real x_lim = x[0];
    real y_lim = x[1];
    real z_lim = x[2];
    real pvap;
    real T_vap = 2835.0;  // Define evaporation temperature
    real L_vap = 3.6e5;  // Define latent heat
    real M_Cu = 0.063546;  // Define molecular weight of copper
    real Pcu = C_P(c, t);  // Get cell pressure
    real T = C_T(c, t);  // Get cell temperature
    real source = 0.0;  // Initialize source term
    real delta_t = CURRENT_TIMESTEP;  // Get current time step size
    pvap = exp(L_vap * (1.0 / T_vap - 1.0 / T) / 56.94);  // Calculate vapor pressure
    real omiga = (pvap * M_Cu) / (pvap * M_Cu + (1.0 - pvap) * 0.14606);  // Calculate mass fraction parameter
    real b = RP_Get_Real("flow-time");  // Get current flow time

    if (b <= 0.000125)
    {
        source = 0.0;  // If time is less than 0.000125, source term is 0
        C_UDMI(c, t, 31) = 0.0;  // Initialize evaporated mass to 0
        C_UDMI(c, t, 32) = 0.0;  // Initialize condensed mass to 0
        C_UDMI(c, t, 33) = 0.0;  // Initialize mass flux term to 0
    }
    else {
        // Check whether in the top-surface ablation region
        if (is_on_top_surface(x_lim, y_lim, z_lim))
        {
            if (T >= T_vap)
            {
                source = 0.06 * 4.0 * 3.1415926 * net(T) / 4.6e6;  // Top-surface evaporation source term
                C_UDMI(c, t, 31) = 0.06 * 4.0 * 3.1415926 * net(T) / 4.6e6;  // Store evaporated mass
                C_UDMI(c, t, 32) = 0.0;  // Condensed mass is 0
                C_UDMI(c, t, 33) = 0.06 * 4.0 * 3.1415926 * net(T) - omiga * Pcu * M_Cu / (8.314 * T * delta_t);
            }
            else
            {
                source = 0.06 * 4.0 * 3.1415926 * net(T) / 4.6e6 - omiga * Pcu * M_Cu / (8.314 * T * delta_t);  // Top-surface net source term
                C_UDMI(c, t, 31) = 0.06 * 4.0 * 3.1415926 * net(T) / 4.6e6;  // Store evaporated mass
                C_UDMI(c, t, 32) = -omiga * Pcu * M_Cu / (8.314 * T * delta_t);  // Store condensed mass
                C_UDMI(c, t, 33) = 0.06 * 4.0 * 3.1415926 * net(T) - omiga * Pcu * M_Cu / (8.314 * T * delta_t);
            }
        }
        // Check whether in the side-surface ablation region
        else if (is_on_side_surface(x_lim, y_lim, z_lim))
        {
            if (T >= T_vap)
            {
                source = 0.06 * 4.0 * 3.1415926 * net(T) / 4.6e6;  // Side-surface evaporation source term
                C_UDMI(c, t, 31) = 0.06 * 4.0 * 3.1415926 * net(T) / 4.6e6;  // Store evaporated mass
                C_UDMI(c, t, 32) = 0.0;  // Condensed mass is 0
                C_UDMI(c, t, 33) = 0.06 * 4.0 * 3.1415926 * net(T) - omiga * Pcu * M_Cu / (8.314 * T * delta_t);
            }
            else
            {
                source = 0.06 * 4.0 * 3.1415926 * net(T) / 4.6e6 - omiga * Pcu * M_Cu / (8.314 * T * delta_t);  // Side-surface net source term
                C_UDMI(c, t, 31) = 0.06 * 4.0 * 3.1415926 * net(T) / 4.6e6;  // Store evaporated mass
                C_UDMI(c, t, 32) = -omiga * Pcu * M_Cu / (8.314 * T * delta_t);  // Store condensed mass
                C_UDMI(c, t, 33) = 0.06 * 4.0 * 3.1415926 * net(T) - omiga * Pcu * M_Cu / (8.314 * T * delta_t);
            }
        }
        else
        {
            source = 0.0;  // If not in the ablation region, source term is 0
            C_UDMI(c, t, 31) = 0.0;  // Evaporated mass is 0
            C_UDMI(c, t, 32) = 0.0;  // Condensed mass is 0
            C_UDMI(c, t, 33) = 0.06 * 4.0 * 3.1415926 * net(T) - omiga * Pcu * M_Cu / (8.314 * T * delta_t);
        }
    }
    C_UDMI(c, t, 30) = source;  // Store total source term value
    dS[eqn] = 0.0;  // Set derivative of source term to 0
    return source;  // Return source term value
}

// Define mass change profile
DEFINE_PROFILE(mass_change, thread, position)
{
    face_t f;
    cell_t c0;
    Thread* t0;
    real b = RP_Get_Real("flow-time");  // Get current flow time

    real x_lim, y_lim, z_lim;
    real T_vap = 2835.0;  // Define evaporation temperature
    real L_vap = 3.6e5;  // Define latent heat
    real M_Cu = 0.063546;  // Define molecular weight of copper
    real delta_t = CURRENT_TIMESTEP;  // Get current time step size

    t0 = THREAD_T0(thread);  // Get thread

    if (b <= 0.000125)
    {
        begin_f_loop(f, thread)
        {
            c0 = F_C0(f, thread);  // Get adjacent cell of face f
            real T = C_T(c0, t0);  // Get cell temperature
            C_UDMI(c0, t0, 33) = 0;  // Initialize a variable to 0
            F_PROFILE(f, thread, position) = 0;  // Set mass change to 0
        }
        end_f_loop(f, thread)
    }
    else
    {
        begin_f_loop(f, thread)
        {
            c0 = F_C0(f, thread);  // Get adjacent cell of face f
            real x[ND_ND];
            C_CENTROID(x, c0, t0);  // Get cell centroid coordinates
            x_lim = x[0];
            y_lim = x[1];
            z_lim = x[2];

            real Pcu = C_P(c0, t0);  // Get cell pressure
            real T = C_T(c0, t0);  // Get cell temperature
            real pvap = exp(L_vap * (1.0 / T_vap - 1.0 / T) / 56.94);  // Calculate vapor pressure
            real omiga = (pvap * M_Cu) / (pvap * M_Cu + (1.0 - pvap) * 0.14606);  // Calculate a parameter

            C_UDMI(c0, t0, 33) = 0.06 * 4.0 * 3.1415926 * net(T);  // Store a variable
            F_PROFILE(f, thread, position) = 0.06 * 4.0 * 3.1415926 * net(T) - omiga * Pcu * M_Cu / (8.314 * T * delta_t);  // Set mass change
        }
        end_f_loop(f, thread)
    }
}
