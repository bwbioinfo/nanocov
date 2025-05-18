// src/plotting/themes/mod.rs
//
// Color themes for plotting coverage data

use plotters::style::RGBColor;

/// Color theme definition for plots
pub struct ColorTheme {
    /// Primary color for main coverage bars
    pub primary: RGBColor,
    /// Secondary/accent color for highlights, mean lines, etc.
    pub accent: RGBColor,
    /// Color for high coverage regions
    pub high: RGBColor,
    /// Color for low coverage regions
    pub low: RGBColor,
    /// Background color
    pub base: RGBColor,
    /// Overlay/box background color 
    pub overlay: RGBColor,
    /// Text color
    pub text: RGBColor,
}

// Catppuccin Latte (light theme)
pub const CATPPUCCIN_LATTE: ColorTheme = ColorTheme {
    primary: RGBColor(30, 102, 245),  // Blue
    accent: RGBColor(136, 57, 239),   // Mauve
    high: RGBColor(234, 83, 83),      // Red
    low: RGBColor(64, 160, 43),       // Green
    base: RGBColor(239, 241, 245),    // Base
    overlay: RGBColor(220, 224, 232), // Overlay
    text: RGBColor(76, 79, 105),      // Text
};

// Catppuccin Frappe (dark theme)
pub const CATPPUCCIN_FRAPPE: ColorTheme = ColorTheme {
    primary: RGBColor(140, 170, 238), // Blue
    accent: RGBColor(186, 187, 241),  // Mauve
    high: RGBColor(231, 130, 132),    // Red
    low: RGBColor(166, 209, 137),     // Green
    base: RGBColor(48, 52, 70),       // Base
    overlay: RGBColor(65, 69, 89),    // Overlay
    text: RGBColor(198, 208, 245),    // Text
};

// Nord Theme
pub const NORD: ColorTheme = ColorTheme {
    primary: RGBColor(94, 129, 172),  // Nord9 (blue)
    accent: RGBColor(180, 142, 173),  // Nord15 (purple)
    high: RGBColor(191, 97, 106),     // Nord11 (red)
    low: RGBColor(163, 190, 140),     // Nord14 (green)
    base: RGBColor(236, 239, 244),    // Nord6 (for light theme)
    overlay: RGBColor(229, 233, 240), // Nord5
    text: RGBColor(46, 52, 64),       // Nord0
};

// Gruvbox Light
pub const GRUVBOX_LIGHT: ColorTheme = ColorTheme {
    primary: RGBColor(69, 133, 136),  // Aqua
    accent: RGBColor(177, 98, 134),   // Purple
    high: RGBColor(204, 36, 29),      // Red
    low: RGBColor(152, 151, 26),      // Green
    base: RGBColor(251, 241, 199),    // Background
    overlay: RGBColor(235, 219, 178), // Light background
    text: RGBColor(60, 56, 54),       // Foreground
};