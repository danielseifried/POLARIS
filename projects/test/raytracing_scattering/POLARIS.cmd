<common>

    <dust_component>    "input/dust_nk/silicate_d03.nk" "plaw" 1.0 3800.0 1e-06 1e-06 -3.5
    <phase_function>    PH_MIE

    <mass_fraction>    0.01

    <nr_threads>    -1

</common>

<task> 1

    <cmd>    CMD_DUST_EMISSION

    <detector_dust_polar nr_pixel = "255*255">    1e-5    1e-3    5    1    0.0    0.0    4.32e+18

    <max_subpixel_lvl>    1

    <path_grid>    "projects/test/raytracing_scattering/grid_3D_sphere_const_T_m1e-5.dat"
    <path_out>    "projects/test/raytracing_scattering/dust/"

</task>

<task> 1

    <cmd>    CMD_DUST_EMISSION

    <detector_dust_polar nr_pixel = "255*255">    1e-5    1e-3    5    1    0.0    0.0    4.32e+18

    <source_star nr_photons = "1e6">    0    0    0    2    4500
    <source_dust nr_photons = "2.1e6">

    <max_subpixel_lvl>    1

    <rt_scattering> 1

    <path_grid>    "projects/test/raytracing_scattering/grid_3D_sphere_const_T_m1e-5.dat"
    <path_out>    "projects/test/raytracing_scattering/dust_rt/"

</task>

<task> 1

    <cmd>    CMD_DUST_SCATTERING

    <detector_dust_mc nr_pixel = "255*255">    1e-5    1e-3    5    0.00    0.00    4.32e+18

    <source_star nr_photons = "1e6">    0    0    0    2    4500
    <source_dust nr_photons = "2.1e6">

    <path_grid>    "projects/test/raytracing_scattering/grid_3D_sphere_const_T_m1e-5.dat"
    <path_out>    "projects/test/raytracing_scattering/dust_mc/"

    <peel_off>    1
    <enfsca>    1

</task>

# this simulation should not be used for gitlab CI/CD as it takes some time
<task> 0

    <cmd>    CMD_DUST_SCATTERING

    <detector_dust_mc nr_pixel = "255*255">    1e-5    1e-3    5    0.00    0.00    4.32e+18

    <source_star nr_photons = "1e10">    0    0    0    2    4500
    <source_dust nr_photons = "2.1e10">

    <path_grid>    "projects/test/raytracing_scattering/grid_3D_sphere_const_T_m1e-5.dat"
    <path_out>    "projects/test/raytracing_scattering/dust_mc_nopeel/"

    <peel_off>    0
    <enfsca>    1

</task>
