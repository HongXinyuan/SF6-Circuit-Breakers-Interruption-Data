// Define geometric parameters (unit: meter)
// Top surface disk parameters
#define TOP_X 0.00234              // X coordinate of top surface center (2.34mm)
#define TOP_RADIUS 0.001           // Top surface radius (1mm)

// Bottom surface disk parameters (only used to define side surface range)
#define BOTTOM_X -0.00766          // X coordinate of bottom surface center (-7.66mm)
#define BOTTOM_RADIUS 0.011        // Bottom surface radius (11mm)

// Side arc parameters (in XY plane)
#define ARC_CENTER_X -0.00766      // X coordinate of arc center (-7.66mm)
#define ARC_CENTER_Y 0.001         // Y coordinate of arc center (1mm)
#define ARC_RADIUS 0.010           // Arc radius (10mm)
#define THETA_MIN 0.0              // Start angle (0 degrees)
#define THETA_MAX 90.0             // End angle (90 degrees)

// Ablation region thickness parameter
#define SURFACE_THICKNESS 0.0001   // Ablation surface thickness (0.1mm)

// Define the region where ablation occurs on the top and side surfaces
#define TOP_SURFACE_MIN_X (TOP_X - SURFACE_THICKNESS/2.0)     // Minimum X of top surface
#define TOP_SURFACE_MAX_X (TOP_X + SURFACE_THICKNESS/2.0)     // Maximum X of top surface

// UDM variable definitions
#define delta_m C_UDMI(c,t,30)     // User-defined memory variable storing copper vapor generation rate at current time step
#define mass_plus C_UDMI(c,t,31)   // Store ablation production term
#define mass_minus C_UDMI(c,t,32)  // Store condensation consumption term

static real cumulative_mass = 0.0;  // Define cumulative mass

// Helper function: check whether a point is in the top-surface ablation region
static int is_on_top_surface(real x, real y, real z)
{
    real r = sqrt(y * y + z * z);  // Distance to X-axis
    int result = 0;

    if (x >= TOP_SURFACE_MIN_X && x <= TOP_SURFACE_MAX_X && r <= TOP_RADIUS)
        result = 1;

    return result;
}

// Helper function: check whether a point is in the side-surface ablation region
static int is_on_side_surface(real x, real y, real z)
{
    // Calculate distance to X-axis
    real r = sqrt(y * y + z * z);

    // Calculate distance from the point to the arc center (in the meridional plane)
    real dx = x - ARC_CENTER_X;
    real dr = r - ARC_CENTER_Y;
    real distance_to_arc_center = sqrt(dx * dx + dr * dr);

    // Check whether the point is near the arc surface (considering thickness)
    int is_near_arc = 0;
    if (fabs(distance_to_arc_center - ARC_RADIUS) <= SURFACE_THICKNESS)
        is_near_arc = 1;

    // Check whether the point is within the side-surface axial range
    // For the side surface, x should be between the bottom and top surfaces
    int is_in_axial_range = 0;
    if (x >= BOTTOM_X && x <= TOP_X)
        is_in_axial_range = 1;

    // Radial range: r should be between arc center Y and bottom radius
    int is_in_radial_range = 0;
    if (r >= ARC_CENTER_Y && r <= BOTTOM_RADIUS)
        is_in_radial_range = 1;

    // Note: since this is a body of revolution, the angular range should be 0-360 degrees,
    // and ablation occurs over the entire side surface.
    // Therefore, no angular restriction is applied to ensure the whole side surface is considered.

    // Return the combined judgment result
    return (is_near_arc && is_in_axial_range && is_in_radial_range);
}

// Helper function: compute the normal direction on the side surface
static void get_side_surface_normal(real x, real y, real z, real* nx, real* ny, real* nz)
{
    // Calculate distance to X-axis
    real r = sqrt(y * y + z * z);

    // Normal direction of the arc in the meridional plane
    real dx = x - ARC_CENTER_X;
    real dr = r - ARC_CENTER_Y;
    real norm = sqrt(dx * dx + dr * dr);

    if (norm > 0.0) {
        // Normal direction (pointing outward from the arc)
        *nx = dx / norm;

        // Components of the radial component in y and z directions (considering 3D)
        if (r > 0.0) {
            *ny = dr * (y / r) / norm;
            *nz = dr * (z / r) / norm;
        }
        else {
            *ny = 0.0;
            *nz = dr / norm;  // When r=0, assume the normal is in the z direction
        }
    }
    else {
        *nx = 1.0;
        *ny = 0.0;
        *nz = 0.0;
    }
}

// Helper function: compute the top-surface normal direction (always along +X axis)
static void get_top_surface_normal(real x, real y, real z, real* nx, real* ny, real* nz)
{
    *nx = 1.0;  // Top-surface normal points along +X axis
    *ny = 0.0;
    *nz = 0.0;
}