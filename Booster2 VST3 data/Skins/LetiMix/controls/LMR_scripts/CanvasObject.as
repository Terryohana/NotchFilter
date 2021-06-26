namespace LM{

// order of parameters used for object coordinates in "o" array 
enum objectCoordsArrayParams
{
  OC_CORNER_X, OC_CORNER_Y, OC_CIRC_CENTER_X, OC_CIRC_CENTER_Y, 
  OC_RADIUS, OC_ARCANGLE_FROM, OC_ARCANGLE_TO, OC_ARCANGLE_FROM_RAD, OC_ARCANGLE_TO_RAD, 
  OC_CORNER_ANGLE_RAD,
  OC_BOT_X, OC_BOT_Y, OC_BOT_Z, OC_TOP_X, OC_TOP_Y, OC_TOP_Z
}

// normal polygons params (polygons for sides)
enum polygonParams{
  POLY_TYPE, POLY_VISIBLE, POLY_COLOR_R, POLY_COLOR_G, POLY_COLOR_B, POLY_COLOR_A, POLY_X1, POLY_Y1, POLY_X2, POLY_Y2, POLY_X3, POLY_Y3, POLY_X4, POLY_Y4, POLY_ANGLENORMALRAD, 
}
// polygons type 0 (POLY_TYPE_OBJECT2D) - simple 2D repeat of object
enum polygon0Params{
  POLY_XOFFSET = POLY_COLOR_A+1, POLY_YOFFSET, POLY_ZOOM, POLY_ROTATION
}
// polygon types
enum polygonTypes{
  POLY_TYPE_OBJECT2D,
  POLY_TYPE_NORMAL,
  POLY_TYPE_DECOR_BEVEL
}

// params in "shadows" array, keeping info about activated shadows
enum shadowParams{
  SH_TYPE, SH_COLOR_R, SH_COLOR_G, SH_COLOR_B, SH_INTENSITY, SH_LIGHT_SOURCE, SH_MAX_BLUR_LAYERS
}

// params in "shadowsPreCalc" array, keeping info about calculated shadows
enum shadowCalcParams{
  SC_FULLY_OVERLAP, SC_COLOR_R, SC_COLOR_G, SC_COLOR_B, SC_OPACITY,
  SC_S1XC, SC_S1YC, SC_S1_RATIO, SC_S2XC, SC_S2YC, SC_S2_RATIO, 
  SC_P1S1X, SC_P1S1Y, SC_P1S2X, SC_P1S2Y, 
  SC_P2S1X, SC_P2S1Y, SC_P2S2X, SC_P2S2Y,
  SC_GFROMX, SC_GFROMY, SC_GTOX, SC_GTOY, SC_GFROMOP, SC_GTOOP, SC_GDISTA, SC_GDISTB, SC_BLURZOOMDELTA, SC_MAXLAYERS
}

// order of params for "scoord" array (used when calculating shadows corner points)
enum shadowCoordsArrayParams{
  SCP_S1_X, SCP_S1_Y, SCP_S2_X, SCP_S2_Y
}

// helper object for dev mode

// generic canvas object for knobs, fader thumbs etc
class CanvasObject {
  int debug_mode = 0;

  bool enable_3d = false;

  // canvas and gui params
  double cw, ch; // original canvas width and height
  double cw2, ch2; // canvas size and half size
  double gui_canvas_x, gui_canvas_y; // left-top coordinates of CANVAS about GUI center
  double gui_x = 0, gui_y = 0; // center coordinates of OBJECT about center of GUI
  double w2, h2; // half of object size and width, ro=radius of roundness
  double xc_offset, yc_offset; // object position relative to canvas center (0,0). Really object is positioned always at (0,0) to enable rotation and scaling, and canvas is shifted by this offset (xc_offset, yc_offset)

  // object parameters
  double xc_orig, yc_orig; // as received - position of center of object relative to top left canvas coords
  double w_orig, h_orig, tallness_orig, elevation_orig, round_orig;
  double w, h, tallness, elevation, round, maxradius; // object properties
  double tallelev; // sum of tallness and elevation
  double elev_percent = 0; // how much object is elevated compared to tallelev

  // rotation params
  double rotation_base_orig, rotation_base, rotation_delta = 0;
  double rotation, rotR; // calculated in degrees and radians
  bool is_rotated = false;

  // tilt params
  double tilt = 0, tiltR = 0, tilt_angle = 0, tilt_angleR; // tilt from -1 to 1, tilt angle from 180 to -180
  double tilt_orig = 0, tilt_angle_orig = 0;

  bool is_perfect_circle = false;
  bool is_flat_on_surface = false;
  
  // array for object shape
  array<array<double>> o(0, {0,0,0,0}); // object coordinates (corner points etc)
  int ncorners_orig = 4, ncorners = 4; // basically equals o.length
  array<array<double>> polygons(); // array for calculated polygons in 3d mode

  // for quicker access to reqular 4-corner objects, we have the first 4 points points in vars
  double o0x,o0y, o1x,o1y, o2x,o2y, o3x,o3y; // coorinates of first four corner points (for quicker access)
  double oc0x,oc0y, oc1x,oc1y, oc2x,oc2y, oc3x,oc3y; // centers of circles of round edges
  double top_x_offset = 0 , top_y_offset = 0, top_zoom = 1; // for 3d mode top surface

  // body colors
  double bodyR, bodyG, bodyB, bodyH, bodyS, bodyL, bodyA; // colors in rgb, hsl and alpha
  double borderR, borderG, borderB, borderA;
  double borderWidth, bodyOpacity, borderOpacity;

  // marker ca
  bool has_marker = false;
  double markerR, markerG, markerB, markerH, markerS, markerL, markerA;
  double markerBorderR, markerBorderG, markerBorderB, markerBorderA = 0;

  // shadows
  array<array<double>> shadows; // to keep params of added shadows
  array<array<double>> shadowsPreCalc; // to keep values of calculated shadows
  bool has_shadows = false; // for quick skipping if no shadows
  array<array<double>> scoord = { {0,0,0,0,0,0,0,0} }; // used while calculation

  // light
  double light_dependency = 0.5; // how much object color (lightness) is affected by light sources
  double top_light_dependency_ratio = 0.1; // for adjusting top surface light dependency (comparing to sides)
  
  // shades and highlights
  double shadeIntensity = 0;
  double highlightIntensity = 0;
  double CL_r, CL_g, CL_b, CL_AngleR, CL_AltR, CL_distance; // combined light source properties
  double ShadeX, ShadeY; // coords for shading gradient
  double ShadeOpStart, ShadeOpEnd, HLOpStart, HLOpEnd; // params for gradient opacity
  bool swap_shade_edges = false, has_shading = false, has_highlight = false;

  // decorated bevel
  double decorBevelHLIntensity, decorBevelShadeIntensity, decorBevelWidth, decorBevelWidth_orig = -1;
  bool has_decor_bevel = false;
  double decorBevelHLPosA, decorBevelHLPosB, decorBevelHLPosAOp, decorBevelHLPosBOp;
  double decorBevelShPosA, decorBevelShPosB, decorBevelShPosAOp, decorBevelShPosBOp;

  //////////////////////////////////////////
  // SETUP AND CALCULATE METHODS
  //////////////////////////////////////////

  // object initialization
  void init(double cw, double ch, double gui_canvas_x, double gui_canvas_y, bool enable_3d = false){
    
    // canvas size and half of it (in pixels)
    this.cw = cw;
    this.ch = ch;
    this.cw2 = cw/2;
    this.ch2 = ch/2;
    this.gui_canvas_x = gui_canvas_x;
    this.gui_canvas_y = gui_canvas_y;
    this.enable_3d = enable_3d;
  }

  // main function to set object shape 
  void setShape(int number_of_angles, double w, double h, double tallness, double elevation, double round, double rot, double tilt = 0, double tilt_angle = 0, bool update_orig = true) {
    // object size, half size and center (in pixels)
    this.w = w;
    this.h = h;
    this.w2 = w/2;
    this.h2 = h/2;
    this.maxradius = this.w2;
    if (this.h2 > this.w2) this.maxradius = this.h2;

    this.tallness = tallness;
    this.elevation = elevation;
    this.tallelev = tallness + elevation;
    if (tallelev != 0) {
      this.elev_percent = elevation/tallelev;  
    } else {
      this.elev_percent = 0;
    }
    
    // roundness (ratio)
    if (round > 1) round = 1;
    if (round < 0) round = 0;
    this.round = round;

    is_perfect_circle = ((w == h) and (round == 1) and (number_of_angles == 4));

    // set object rotation
    setRotation(rot, false);

    // set number of angles
    if (number_of_angles < 3) number_of_angles = 3;
    if (number_of_angles > 100) number_of_angles = 100;
    this.ncorners = number_of_angles;

    // set tilt
    this.tilt = tilt;
    this.tilt_angle = tilt_angle;
    tiltR = tilt*pi/180;
    tilt_angleR = tilt_angle*pi/180;

    // calculate objects corner points and keep in "o" array
    calcCornerPoints(o);

    // remember initial settings to be able to change them
    if (update_orig) {
      this.w_orig = w;
      this.h_orig = h;
      this.tallness_orig = tallness;
      this.round_orig = round;
      this.rotation_base_orig = rotation;
      this.elevation_orig = elevation;
      this.ncorners_orig = number_of_angles;
      this.tilt_orig = tilt;
      this.tilt_angle_orig = tilt_angle;
    }

    is_flat_on_surface = (tallelev == 0);
  }

  // set general
  void setRotation(double rot, bool redraw = true){
    // object rotation
    this.rotation_base = rot;

    rotation = rotation_base + rotation_delta;
    rotR = (rotation)*pi/180; // in Radians
    if (rotation == 0) is_rotated = false; else is_rotated = true;
    if (redraw) {
      calcShadows();
      calcShadingAndBevel();
      calcPolygons();
      
    }
  }

  // usable for knobs and other params rotating relatively to rotation
  void setRotationDelta(double rotation_delta, bool redraw = true){
    this.rotation_delta = rotation_delta;
    setRotation(rotation_base, redraw);
  }

  // on render settings change
  void updateRenderSettings() { 
    calcShadows();
    calcShadingAndBevel();
    calcPolygons();
  }

  // calculate objects corner points, radiuses etc.
  void calcCornerPoints(array<array<double>> &out coords){

    // put all corner points into array
    coords.resize(ncorners);
    scoord.resize(ncorners); // for shadow data

    // calculate coords for all corner points
    if (ncorners == 4) {
        // for simple traditional objects with 4 angles
        // roundness of edges (in pixels)
        double ro;
        if (w > h) ro = round * h/2.0; else ro = round * w/2.0;

        coords.resize(ncorners+1);
        coords[0] = {-w2, -h2, -w2+ro, -h2+ro, ro,  -180, -90,  -pi, -pi2,  0};
        coords[1] = {w2, -h2, w2-ro, -h2+ro, ro,    -90, 0,     -pi2, 0,    0};
        coords[2] = {w2, h2, w2-ro, h2-ro, ro,      0, 90,      0, pi2,     0};
        coords[3] = {-w2, h2, -w2+ro, h2-ro, ro,    90, 180,    pi2, pi,    0};
        
      } else {
        // when object has arbitrary number of angles
        // generate auto-coordinates for this figure

        // determine starting angle
        //double start_angle_offsetR = 0; // start from top
        double cur_angleR = -pi2; //+ start_angle_offsetR;

        // determine for which angle we rotate each step
        double angle_deltaR = (360.0f/double(ncorners))*pi/180;
        double angle_deltaR2 = angle_deltaR/(2);

        // calculate max distance till point in angle
        double distance = w2;
        // if (ncorners == 4) distance *= sqrt(2); // for "normal" figure
        double hwratio = h/w; // how height relates to width

        for(int n=0;n<ncorners;n++) {
          double co = cos(cur_angleR), si = sin(cur_angleR);
          
          // corner points
          double ex = co*distance, ey = si*distance*hwratio;

          // rounding centers and radiuses     
          double radius = (distance*round)/(pow(hwratio, 3));
          double radius_delta_x = distance - radius;
          double radius_delta_y = distance*hwratio - radius;
          double ecx = co*radius_delta_x;
          double ecy = si*radius_delta_y;

          // start and end angles of circle
          double angle_startR = (cur_angleR-angle_deltaR2);
          double angle_endR = (cur_angleR+angle_deltaR2);

          // add corner point
          // if (lm_cur_angleR > pi) lm_cur_angleR -= twopi;
          coords[n] = {ex, ey, ecx, ecy, radius, angle_startR*180/pi, angle_endR*180/pi, angle_startR, angle_endR, 0};

          // increase angle to next point
          cur_angleR += angle_deltaR;
        }
      }

      // calculate angles to corners
      for(int n=0;n<ncorners;n++) {
        // double ocx = coords[n][OC_CIRC_CENTER_X], ocy = coords[n][OC_CIRC_CENTER_Y];
        double ox = coords[n][OC_CORNER_X], oy = coords[n][OC_CORNER_Y];
        double angleR = de(atan2(ox, -oy));
        coords[n][OC_CORNER_ANGLE_RAD] = angleR;
      }
    
      // add tilt
      if (tiltR != 0)
      for(int n=0;n<ncorners;n++) {
        double ox = coords[n][OC_CORNER_X], oy = coords[n][OC_CORNER_Y];
        double cAngleR = coords[n][OC_CORNER_ANGLE_RAD]; 
        // find triangle angle and sides
        double angle_diffR_orig = angDiffR(tilt_angleR, cAngleR);
        double c_orig = sqrt((ox*ox) + (oy*oy)); // hypothenuse
        double b_orig = de(c_orig*cos(angle_diffR_orig)); // cathet b
        double a_orig = de(c_orig*sin(angle_diffR_orig)); // cathet a
        
        // now we apply tilting
        double b = de(cos(tiltR)*b_orig);
        double angle_diffR = atan(a_orig/b); // angle difference betwee tilt angle and point angle
        double angleR = tilt_angleR - angle_diffR;
        double c = b/cos(angle_diffR);
        if (b==0) c = abs(a_orig);

        // string text = "";
        // text += "c_orig: " + c_orig + " b_orig: " + b_orig + " a_orig:" + a_orig + " b: " + b;
        // status(text);

        coords[n][OC_CORNER_X] = sin(angleR)*c;
        coords[n][OC_CORNER_Y] = -cos(angleR)*c;

        coords[n][OC_CIRC_CENTER_X] = coords[n][OC_CORNER_X];
        coords[n][OC_CIRC_CENTER_Y] = coords[n][OC_CORNER_Y];
        coords[n][OC_RADIUS] = 0;
        // TODO: tilting of rounded corners?

        // double a = tan(angle)*b;
        //status("angle: " + angle_diffR_orig*180/pi + " angleNew: " + angle_diffR*180/pi + " angleR: " + angleR*180/pi + " c_orig: " + c_orig + " b_orig: " + b_orig + " a_orig:" + a_orig + " b: " + b);
      }

    if (enable_3d) {
      // for 3d objects add 3d info


    }

    // for simple objects
    if (ncorners == 4) {
      // for quicker access copy first 4 array values to vars
      o0x = coords[0][OC_CORNER_X]; o0y = coords[0][OC_CORNER_Y]; 
      oc0x = coords[0][OC_CIRC_CENTER_X];  oc0y = coords[0][OC_CIRC_CENTER_Y];
      o1x = coords[1][OC_CORNER_X]; o1y = coords[1][OC_CORNER_Y]; 
      oc1x = coords[1][OC_CIRC_CENTER_X];;  oc1y = coords[1][OC_CIRC_CENTER_Y];
      o2x = coords[2][OC_CORNER_X]; o2y = coords[2][OC_CORNER_Y];
      oc2x = coords[2][OC_CIRC_CENTER_X];  oc2y = coords[2][OC_CIRC_CENTER_Y];
      o3x = coords[3][OC_CORNER_X]; o3y = coords[3][OC_CORNER_Y]; 
      oc3x = coords[3][OC_CIRC_CENTER_X];  oc3y = coords[3][OC_CIRC_CENTER_Y];
    }

  }

  void setBorder(string borderColor, double borderOpacity, double borderWidth){
    // border color and width
    borderR = double(parseInt(borderColor.substr(1,2),16))/255.0;
    borderG = double(parseInt(borderColor.substr(3,2),16))/255.0;
    borderB = double(parseInt(borderColor.substr(5,2),16))/255.0;
    borderA = borderOpacity;
    this.borderWidth = borderWidth;
  }

  void setBodyColor(string bodyColor, double bodyOpacity){
    // body color
    bodyR = double(parseInt(bodyColor.substr(1,2),16))/255.0;
    bodyG = double(parseInt(bodyColor.substr(3,2),16))/255.0;
    bodyB = double(parseInt(bodyColor.substr(5,2),16))/255.0;
    bodyA = bodyOpacity;
    convertRGBtoHSL(bodyR, bodyG, bodyB, bodyH, bodyS, bodyL);
  }

  void setBodyColorHSL(double h, double s, double l){
    bodyH = h; bodyS = s; bodyL = l;
    convertHSLtoRGB(bodyH, bodyS, bodyL, bodyR, bodyG, bodyB);
  }

  // add marker information
  void setMarkerColors(string markerColor, double markerOpacity, string markerBorderColor, double markerBorderOpacity) {
    
    has_marker = true;
    
    // prepare colors
    markerR = double(parseInt(markerColor.substr(1,2),16))/255.0;
    markerG = double(parseInt(markerColor.substr(3,2),16))/255.0;
    markerB = double(parseInt(markerColor.substr(5,2),16))/255.0;
    convertRGBtoHSL(markerR, markerG, markerB, markerH, markerS, markerL);
    markerA = markerOpacity;

    markerBorderR = double(parseInt(markerBorderColor.substr(1,2),16))/255.0;
    markerBorderG = double(parseInt(markerBorderColor.substr(3,2),16))/255.0;
    markerBorderB = double(parseInt(markerBorderColor.substr(5,2),16))/255.0;
    markerBorderA = markerBorderOpacity;
  }


  // add shading to emulate light
  void setShadeAndHighlight(double shade_intensity = 0, double highlight_intensity = 0){
    shadeIntensity = shade_intensity;
    highlightIntensity = highlight_intensity;
  }

  // set decorated bevel (mostly for knobs)
  void setDecorBevel(double width = 0, double hl = 0, double shade = 0) {
    this.decorBevelHLIntensity = hl;
    this.decorBevelShadeIntensity = shade;
    this.decorBevelWidth = width;
    this.has_decor_bevel = (width>0) and ((hl>0) or (shade>0));
    if (decorBevelWidth_orig == -1) {
      decorBevelWidth_orig = decorBevelWidth;
    }
  }


  // set object center position (in left-top coordinates of original canvas)
  void setCenterPos(double xc_orig, double yc_orig){
    this.xc_orig = xc_orig;
    this.yc_orig = yc_orig;
    this.xc_offset = xc_orig-cw2; // we offset canvas keeping object at 0,0
    this.yc_offset = yc_orig-ch2;
    gui_x = gui_canvas_x + xc_orig; // calculate center coordinates of object about center of gui
    gui_y = gui_canvas_y + yc_orig;
    calcShadows();
    calcShadingAndBevel();
    calcPolygons();
  }

  void addShadow(int type, string color, double opacity, int light_source, int max_blur_layers = 10){
    int n = shadows.length;
    shadows.resize(n+1);
    shadowsPreCalc.resize(n+1);

    // add shadow params to shadows array
    double r, g, b;
    r = double(parseInt(color.substr(1,2),16))/255.0;
    g = double(parseInt(color.substr(3,2),16))/255.0;
    b = double(parseInt(color.substr(5,2),16))/255.0;

    shadows[n] = { type + 0.0f, r, g, b, opacity, light_source + 0.0f, max_blur_layers + 0.0f};
    has_shadows = true; 
  }

  // calculate light properties (angles, intensity) for current object position
  void calcLight(double x, double y, double z, array<double>@ ls, double &out LangleR, double &out LaltR, double &out Ldistance, double &out Lsize, double &out Lintensity, double &out Lx, double &out Ly, double &out Lz){

    // calculate new light position (relative to current object)
    Lx = ls[LS_X] - x; Ly = ls[LS_Y] - y; Lz = ls[LS_Z] - z;

    // calculate current light angle in XY space
    LangleR = atan2(Lx,-Ly);

    // vector length to light source
    Ldistance = sqrt(Lx*Lx + Ly*Ly + Lz*Lz); 

    // calculate new angle in Z space
    LaltR = atan2(Lz, sqrt(Ldistance*Ldistance - Lz*Lz));

    // in future we can also modify other light params depending on light distance
    Lsize = ls[LS_SIZE];

    Lintensity = ls[LS_INTENSITY]*(1/pow(Ldistance/Light.ref_intensity_distance, 2));

    //Lintensity = 10*ls[LS_INTENSITY]*(1/pow(Ldistance/100, 2));

    // calculate distance between object and light source
    // ls[LS_DISTANCE]
  }

  // calculate shadows for current object position
  void calcShadows(){
    // clear precalculated shadows data
    shadowsPreCalc.resize(0);

    if ((!has_shadows) or (is_flat_on_surface)) return;

    // vars to keep current light info
    double LangleR, LaltR, Ldistance, Lsize, Lintensity, Lx, Ly, Lz, LsizeRatio;

    // get number of shadows (one "shadow" can support multiple light sources)
    int nshadows = int(shadows.length); 
    if (nshadows == 0) return;

    // get number of light sources
    array<array<double>>@ lights = Light.sources;
    int nlights = int(lights.length);
    if (nlights == 0) return;

    // get object max size
    double objmaxcenterdistance = w2; if (h2>w2) objmaxcenterdistance = h2;

    
    // parse through all "shadows" and precalculate them
    for (int sh_n = 0; sh_n < int(shadows.length); sh_n++) {

      array<double>@ sh = shadows[sh_n];
      double orig_sh_intensity = sh[SH_INTENSITY];
      if (orig_sh_intensity == 0) continue;

      int sh_light_source = rint(sh[SH_LIGHT_SOURCE]);
      int ls_start = 0, ls_end = nlights-1;
      if (sh_light_source >=0) { ls_start = sh_light_source; ls_end = sh_light_source; }

      // step though all (or only required) light sources
      for(int ls_n = ls_start; ls_n <= ls_end; ls_n++) {
        array<double>@ ls = lights[ls_n];
        if (ls[LS_ENABLED] == 0) continue; // skip if light source is off
        //echo("Light source "+ls_n);

        // calculate LangleR, LaltR, Ldistance, Lsize, Lintensity, Lx, Ly, Lz for current light source and object position
        // we get these angles supposing like object is centered at 0,0 point
        calcLight(gui_x, gui_y, 0, ls, LangleR, LaltR, Ldistance, Lsize, Lintensity, Lx, Ly, Lz);

        //echo("Light source: "+m + " LangleR: " + LangleR + " LaltR: " + LaltR + " Ldistance: " + Ldistance);

        // for each shadow object we calculate three things:
        // 1) shadow of the bottom 2) shadow of top of object
        // 3) connection points between them

        // calculate SHADOW 1 (bottom of object)

        double distance_till_light = (Ldistance - elevation); if (distance_till_light < 1) distance_till_light = 1;

        // calculate s1_ratio aware of light size
        // double tan_a = (w*0.5 - Lsize*0.5)/distance_till_light; // we suppose light is on top
        // double half_shadow_increase = tan_a*elevation; // length of added shadow area on one side
        // double s1_ratio_ls = (2*half_shadow_increase + w)/w;
        // if (s1_ratio_ls < 0.00001) s1_ratio_ls = 0.00001; // can be removed some day

        // calculate s1_ratio ignoring light size
        double s1_ratio = (2*((w*0.5)/distance_till_light)*elevation + w)/w;

        double s1w = w*s1_ratio, s1h = h*s1_ratio; // shadow 1 body size

        // distance till shadow 1 depends on elevation of object
        double s1_distance = elevation/tan(LaltR);

        // calculate shadow 1 offset
        double LAsin = sin(LangleR - rotR), LAcos = cos(LangleR - rotR);
        double s1xc = -LAsin*s1_distance, s1yc = LAcos*s1_distance; // center of shadow

        // now calculate SHADOW 2 (top of object)
        distance_till_light = (Ldistance - tallelev); if (distance_till_light < 1) distance_till_light = 1;

        // calculate s2_ratio_ls beeing aware of light size
        // tan_a = (w*0.5 - Lsize*0.5)/distance_till_light; // we suppose light is on top
        // half_shadow_increase = tan_a*tallelev; // length of added shadow area on one side
        // double s2_ratio_ls = (2*half_shadow_increase + w)/w;

        // calculate s2_ratio ignoring light size
        double s2_ratio = (2*((w*0.5)/distance_till_light)*tallelev + w)/w;
        

        if (s2_ratio < 0.00001) s2_ratio = 0.00001; // can be removed some day
        if (s2_ratio < (s1_ratio+0.1)) s2_ratio = s1_ratio+0.1; // can be removed some day

        double s2w = w*s2_ratio, s2h = h*s2_ratio; // calculate shadow1 body size

        // distance till shadow depends on tallness of object and light altitude (vertical angle)
        double s2_distance = tallelev/tan(LaltR);   
        double s2xc = -LAsin*s2_distance, s2yc = LAcos*s2_distance; // center of shadow

        // now we know two shadows (top and bottom)
        double scAngR = atan2(s2xc-s1xc,s2yc-s1yc); // angle between shadow centres
        double scAngRC = (scAngR-pi2)*-1; if (scAngRC<0) scAngRC += 2*pi; // convert shadow angle to angles used in Arc

        //string text = "";

        // now step though all corner points and find connection points between two shadows
        // we also find the far-est point of shadow
        double maxAdiffA=0, maxAdiffB=0;
        int pAn = 1, pBn = 1;
        for (int n=0;n<ncorners;n++){
          array<double>@ r = o[n];
          double ex = r[OC_CORNER_X], ey = r[OC_CORNER_Y];
          double ecx = r[OC_CIRC_CENTER_X], ecy = r[OC_CIRC_CENTER_Y], ro = r[OC_RADIUS];
          double afromR = r[OC_ARCANGLE_FROM]*pi/180, atoR = r[OC_ARCANGLE_TO]*pi/180, aRangeR = (abs(atoR-afromR)), aMidR = (afromR+(atoR-afromR)/2);
          double s1ccx, s1ccy, s2ccx, s2ccy;

          s1ccx = ecx + ro*cos(aMidR);
          s1ccy = ecy + ro*sin(aMidR);

          s2ccx = s1ccx*s2_ratio + s2xc;
          s2ccy = s1ccy*s2_ratio + s2yc;
          s1ccx = s1ccx*s1_ratio + s1xc;
          s1ccy = s1ccy*s1_ratio + s1yc;

          // calculate angle between points
          double angR = atan2((s2ccx-s1ccx),(s2ccy-s1ccy));      
          double angRC = (angR-pi2)*-1; if (angRC<0) angRC += 2*pi; // convert to angles used in Arc

          // add 90deg shift
          double aR = angRC+pi2; if (aR > twopi) aR -= twopi;

          // flip if on the other side
          double maxAllowedAngleR = aRangeR; 
          // if (n == 5) {
          //   text += "aR: "+(aR*180/pi)+", aMidR: " +(aMidR*180/pi);
          //   text += ", angDiff: " + (angDiffR(aR, aMidR)*180/pi);
          // }
          if (abs(angDiffR(aR, aMidR)) > maxAllowedAngleR) aR -= pi;

          // limit to points inside rounding angle
          double cmp;
          cmp = angDiffR(aR, afromR);
          if (((cmp < 0) and (cmp>-maxAllowedAngleR)) or (cmp > maxAllowedAngleR)) aR = afromR;
          cmp = angDiffR(aR, atoR);
          if ((cmp > 0) and (cmp<maxAllowedAngleR)) aR = atoR;

          // ctx.path.Clear();
          // ctx.path.MoveTo(ecx*s2_ratio + s2xc , ecy*s2_ratio + s2yc);
          // ctx.path.Arc(ecx*s2_ratio + s2xc , ecy*s2_ratio + s2yc, ro*s2_ratio, (aR*180/pi), (aR*180/pi));
          // ctx.StrokePath();

          // correct point of round angle
          s1ccx = ecx + ro*cos(aR);
          s1ccy = ecy + ro*sin(aR);

          s2ccx = s1ccx*s2_ratio + s2xc;
          s2ccy = s1ccy*s2_ratio + s2yc;
          s1ccx = s1ccx*s1_ratio + s1xc;
          s1ccy = s1ccy*s1_ratio + s1yc;

          // re-calculate angle between points
          angR = atan2((s2ccx-s1ccx),(s2ccy-s1ccy));
          double angDiff = scAngR - angR;
          if (angDiff>pi) angDiff-=2*pi;
          if (angDiff<-pi) angDiff+=2*pi;

          // determine maximum angle difference
          if (angDiff < maxAdiffA) {
            maxAdiffA = angDiff;
            pAn = n;
          }
          if (angDiff > maxAdiffB) {
            maxAdiffB = angDiff;
            pBn = n;
          }

          scoord[n] = {s1ccx, s1ccy, s2ccx, s2ccy, 1};
          // text += " ang["+n+"]: "+(angDiff*180/pi); 
        }

        // now we now corners that we should connect
        // scoord[pAn][4] = 3;
        // scoord[pBn][4] = 3;
        if (s2_ratio < s1_ratio) {
          int t = pAn;
          pAn = pBn;
          pBn = t;
        }

        // calc angles to detect if shadow connection is pointing "inside" shadow 2 (to prevent path clipping)
        bool shadows_fully_overlap = false;
        {
          bool shadow1_probably_inside = false;
          double xgrow = s2w-s1w, ygrow = s2h-s1h;
          double xoffs = s2xc-s1xc, yoffs = s2yc-s1yc;
          // text += "xgrow: " + xgrow + ", xoffs: " + xoffs + ", ygrow: " + xgrow + ", yoffs: " + xoffs;
          if ((abs(xoffs)<=(xgrow/2)) and (abs(yoffs)<=(ygrow/2))) shadow1_probably_inside = true;

          // if basic calculation show that shadow 1 is probably inside shadow 2, do more calculations
          if (shadow1_probably_inside) {

            double a, b, c, ang1, ang2, p0, p1, p2, p3, div;

            p0 = scoord[pAn][2]; p1 = scoord[pAn][3]; p2 = scoord[pBn][2]; p3 = scoord[pBn][3];

            // calculate angle of triangle between shadow 2 center and shadow 2 edge points
            a = sqrt(pow(s2xc - p0, 2) + pow(s2yc - p1, 2));
            b = sqrt(pow(p2 - p0, 2) + pow(p3 - p1, 2));
            c = sqrt(pow(s2xc - p2, 2) + pow(s2yc - p3, 2));
            div = (2*a*b); // if (div == 0) div = 1;
            if (abs((a+c)-b) < 0.005)  {
              ang1 = 0; // flat triangle
            } else {
              ang1 = acos((a*a + b*b - c*c)/div);
            }
            
            // do the same for triangle between shadow 2 center, one of edges and shadow 1 point
            p0 = scoord[pAn][2]; p1 = scoord[pAn][3]; p2 = scoord[pAn][0]; p3 = scoord[pAn][1];

            a = sqrt(pow(s2xc - p0, 2) + pow(s2yc - p1, 2));
            b = sqrt(pow(p2 - p0, 2) + pow(p3 - p1, 2));
            c = sqrt(pow(s2xc - p2, 2) + pow(s2yc - p3, 2));
            div = (2*a*b); //if (div == 0) text += " div 0a ";
            ang2 = acos((a*a + b*b - c*c)/div);
            //text +=" ang2: "+(ang2*180/pi);

            if (ang2 < ang1) {
              shadows_fully_overlap = true;
            }
          }
        }

        // status_custom_text = " LangleR: " + (LangleR*180/pi) + " LaltR: " + (LaltR*180/pi) + " s2xc: "+s2xc + " s2yc: "+s2yc + " s2_ratio: "+s2_ratio;

        double fdraw_only_top_shadow = 0; if (shadows_fully_overlap) fdraw_only_top_shadow = 1;
        

        // calculate gradient line to draw intensity
        double d = objmaxcenterdistance; // how big the object is (max w or h)

        // gradient start and end points
        double co = cos(LangleR+rotR), si = sin(LangleR+rotR);
        double g_from_x = s1xc - d*si*s1_ratio, g_from_y = s1yc - d*co*s1_ratio, g_to_x = s2xc + d*si*s2_ratio, g_to_y = s2yc + d*cos(LangleR+rotR)*s2_ratio;
        
        // distance from start and end of gradient to light source
        double distA = Ldistance + sqrt(s1xc*s1xc + s1yc*s1yc) - d*s1_ratio;
        double distB = Ldistance + sqrt(s2xc*s2xc + s2yc*s2yc) + d*s2_ratio;
        // set intensity for start and end of gradient
        
        double orig_light_intensity = ls[LS_INTENSITY];
        if (tallelev < 10) orig_light_intensity *= tallelev/10;
        double g_from_op = orig_light_intensity * pow((1/pow(distA/Light.ref_intensity_distance, 2)), 1.2);
        double g_to_op = orig_light_intensity * pow((1/pow(distB/Light.ref_intensity_distance, 2)), 1.2);
        


        //status("orig_sh_intensity: "+orig_sh_intensity + " g_from_op: " + g_from_op + " g_to_op: " + g_to_op);

        // determine how much zoom we apply for "blur" effect
        // double objSize = sqrt(pow((scoord[pAn][SCP_S1_X] - scoord[pBn][SCP_S1_X]),2) + pow((scoord[pAn][SCP_S1_Y] - scoord[pBn][SCP_S1_Y]),2))/s1_ratio;

        double objSize = w; // if (h>w) objSize = h;
        LsizeRatio = 100*(Lsize/2 - objSize/2)/Ldistance;
        //status("Lsize: " + Lsize + " objSize: " + objSize + " LsizeRatio: "+LsizeRatio);

        double blur_intensity = floor(LsizeRatio);
        if (blur_intensity < 0) blur_intensity = 0;
        if (blur_intensity > 10) blur_intensity = 10;
        if (tallelev < 10) blur_intensity *= tallelev/10; // fix to remove blur for small elements

        // if (ls_n == 0) status("LsizeRatio: " + LsizeRatio + " blur_intensity: "+blur_intensity);

        double blurMaxLayers = floor(blur_intensity+1);
        if (blurMaxLayers > sh[SH_MAX_BLUR_LAYERS]) blurMaxLayers = sh[SH_MAX_BLUR_LAYERS];
        double blurZoomDelta = (blurMaxLayers+1)*0.005;
        
        // TODO LIMIT blurMaxLayers from KUIML

        // status("blurZoomDelta: " + blurZoomDelta + " blurMaxLayers: " + blurMaxLayers);

        // insert calculated shadow params into array
        shadowsPreCalc.insertLast({ fdraw_only_top_shadow, sh[SH_COLOR_R], sh[SH_COLOR_G], sh[SH_COLOR_B], orig_sh_intensity, s1xc, s1yc, s1_ratio, s2xc, s2yc, s2_ratio, scoord[pAn][SCP_S1_X], scoord[pAn][SCP_S1_Y], scoord[pAn][SCP_S2_X], scoord[pAn][SCP_S2_Y], scoord[pBn][SCP_S1_X], scoord[pBn][SCP_S1_Y], scoord[pBn][SCP_S2_X], scoord[pBn][SCP_S2_Y], g_from_x, g_from_y, g_to_x, g_to_y, g_from_op, g_to_op, distA, distB, blurZoomDelta, blurMaxLayers });

      } // end of cycling through light sources
    } // end of cycling through shadows


    // echo("shadows: "+shadowsPreCalc.length);
    // echo_flush(lm_text3);

  } // eng of calcShadows method


  // calculate shading positions
  void calcShadingAndBevel(){
    double CL_x = 0, CL_y = 0, CL_z = 0; // vector to combined light source
    CL_r = 0; CL_g = 0; CL_b = 0;

    // pass through active light sources
    int nlights = int(Light.sources.length);
    int nlights_active = 0;
    
    double total_intensity = 0;
    for(int ls_n = 0; ls_n < nlights; ls_n++) {
      
      array<double>@ ls = Light.sources[ls_n];
      if (ls[LS_ENABLED] == 0) continue; // skip if light source is off
      
      nlights_active++;

      double LangleR, LaltR, Ldistance, Lsize, Lintensity, Lx, Ly, Lz;

      // calculate angle to the light source for the current object
      calcLight(gui_x, gui_y, tallelev, ls, LangleR, LaltR, Ldistance, Lsize, Lintensity, Lx, Ly, Lz);

      CL_x += Lx*Lintensity;
      CL_y += Ly*Lintensity;
      CL_z += Lz*Lintensity;

      CL_r += ls[LS_COLOR_R]*Lintensity;
      CL_g += ls[LS_COLOR_G]*Lintensity;
      CL_b += ls[LS_COLOR_B]*Lintensity;

      total_intensity += Lintensity;
    }
    CL_r /= total_intensity;
    CL_g /= total_intensity;
    CL_b /= total_intensity;

    // coordinates of combined light source
    CL_x /= total_intensity;
    CL_y /= total_intensity;
    CL_z /= total_intensity;

    // combined light source properties
    CL_AngleR = atan2(CL_x, -CL_y);
    // vector length to light source
    CL_distance = sqrt(CL_x*CL_x + CL_y*CL_y + CL_z*CL_z); 
    CL_AltR = atan2(CL_z, sqrt(CL_distance*CL_distance - CL_z*CL_z));  

    // get point coords on edges closer to light
    ShadeX = maxradius * sin(CL_AngleR);
    ShadeY = -maxradius * cos(CL_AngleR);

    //status("" + (CL_AngleR*180/pi) + ", " + (CL_AltR*180/pi) + ", " + ShadeX + ", " + ShadeY);

    // intensity of light in center of gui on tallelev of object
    double LI_center = calcPolygonLightIntensity(0, 0, tallelev, 0, 90); 
    // status(""+LI_center);

    // status("CL "+CL_r+" "+CL_g+" "+CL_b+" ");
    double li_a = calcPolygonLightIntensity(gui_x + ShadeX, gui_y + ShadeY, tallelev, 0, 90);
    double li_b = calcPolygonLightIntensity(gui_x - ShadeX, gui_y - ShadeY, tallelev, 0, 90);
    
    ShadeOpStart = 0;
    ShadeOpEnd = 0;
    if (shadeIntensity > 0) {
      // calc shade
      double shade_mod = shadeIntensity;
      double shade_max = 0.7;

      double LI_diff_a = li_a - LI_center; // how much lighter is closer side then light in center
      double LI_diff_b = LI_center - li_b; // how much darker is further side then light in center
      
      if (LI_diff_a < 0) ShadeOpStart = LI_diff_a*shade_mod/3;
      ShadeOpStart /= LI_center;
      if (ShadeOpStart > shade_max) ShadeOpStart = shade_max;
      
      if (LI_diff_b > 0) { ShadeOpEnd = LI_diff_b*shade_mod; }
      ShadeOpEnd /= LI_center;
      if (ShadeOpEnd > shade_max) ShadeOpEnd = shade_max;
    }

    // calc highlight
    HLOpStart = 0;
    HLOpEnd = 0;
    if (highlightIntensity > 0) {
      double highlight_mod = 0.1 * highlightIntensity;
      double HL_intensity_mod = pow(total_intensity, 0.1);
      if (HL_intensity_mod < 1) HL_intensity_mod = 1;

      double highlight_end_mod = 5*highlight_mod;
      double highlight_max = 0.8;
      double LI_center_half = LI_center*0.5;

      double HLI_diff_a = li_a - LI_center_half; // how much lighter is closer side then light in center
      double HLI_diff_b = li_a - li_b; // how much lighter is further side then back side

      if (HLI_diff_a > 0) { HLOpStart = HLI_diff_a*highlight_mod*HL_intensity_mod; }
      HLOpStart /= LI_center_half;
      if (HLOpStart > highlight_max) HLOpStart = highlight_max;
      
      if (HLI_diff_b < 0) {
          HLOpEnd = HLOpStart;
        } else {
          HLOpEnd = HLOpStart - (HLI_diff_b/li_a)*highlight_end_mod;
        }
      if (HLOpEnd > highlight_max) HLOpEnd = highlight_max;

      // if light intensity in center is bigger than on edges
      // if ((LI_center > li_a) and (LI_center > li_b)) {
      //   RadialShade = true;
      //   HLOpStart = highlight_max;
      //   HLOpEnd = 0;
      // } else {
      //   RadialShade = false;
      // }
    }

    has_shading = ((ShadeOpStart + ShadeOpEnd) > 0);
    has_highlight = ((HLOpStart + HLOpEnd) > 0);

    // calc decorated bevel
    // this is all slightly rude
    // and maybe should be done better some day
    if (has_decor_bevel) {
      double anglenval = pow(sin(CL_AltR),2); // changes from 0 to 1 depending on altitude of light source (1 is 90deg on top)

      // calc highlighted part of bevel
      double hlmod = decorBevelHLIntensity;
      double maxval = pow(1-(anglenval/1.1), 1.3);  
      decorBevelHLPosA = 0;   decorBevelHLPosAOp = hlmod*maxval;
      decorBevelHLPosB = 0.5;   
      // decorBevelHLPosBOp = hlmod*maxval*pow(anglenval,1.3); 
      decorBevelHLPosBOp = 0;

      // calc shaded part of bevel
      double shmod = decorBevelShadeIntensity*0.5;
      double minval = anglenval/9;
      decorBevelShPosA = 0.5;   decorBevelShPosAOp = shmod*minval;
      decorBevelShPosB = 1;   decorBevelShPosBOp = shmod*(minval+maxval);
    }
  } // end of calcShadingAndBevel method

  void calcPolygons(){
    if (!enable_3d) return;

    // get number of light sources
    // array<array<double>>@ lights = Light.sources;
    // int nlights = int(lights.length);

    double cr, cg, cb, ca = bodyA; // to keep colors
    double visibility = 1;  // polygon visibility
    double LangleR, LaltR, Ldistance, Lsize, Lintensity, Lx, Ly, Lz; // for light calculation

    // reset polygons array
    polygons.resize(0);

    // get camera properties
    double cam_persp = Camera.perspective;
    double cam_dist = Camera.distance;
    double cam_x_offs = Camera.x_offset;
    double cam_y_offs = Camera.y_offset;

    // CALCULATE BODY BOTTOM
    double from_surface = elevation;
    // surface gets bigger closer to camera
    double objsizechange = (1/cam_dist)*(cam_dist - from_surface);
    if (objsizechange<0.00001) objsizechange = 0.00001;
    double bot_zoom = 1/objsizechange;
    double bot_x_offset = gui_x*cam_persp*0.001*from_surface*0.01*bot_zoom + cam_x_offs*from_surface*0.01;
    double bot_y_offset = gui_y*cam_persp*0.001*from_surface*0.01*bot_zoom + cam_y_offs*from_surface*0.01;


    // add bottom polygon
    visibility = 1;
    convertHSLtoRGB(bodyH, bodyS, bodyL*0.5, cr,cg,cb);
    // polygons.insertLast({POLY_TYPE_OBJECT2D, visibility, cr, cg, cb, ca, bot_x_offset, bot_y_offset, bot_zoom, rotation});

    // DRAW BOTTOM BORDER
    // ctx.settings.set_lineWidth(borderWidth);
    // ctx.source.SetRGBA (borderR,borderG,borderB, borderA);
    // ctx.StrokePath();

    // CALCULATE BODY TOP
    from_surface = tallelev;
    objsizechange = (1/cam_dist)*(cam_dist - from_surface);
    if (objsizechange<0.00001) objsizechange = 0.00001;
    top_zoom = 1/objsizechange;
    top_x_offset = gui_x*cam_persp*0.001*from_surface*0.01*top_zoom + cam_x_offs*from_surface*0.01;
    top_y_offset = gui_y*cam_persp*0.001*from_surface*0.01*top_zoom + cam_y_offs*from_surface*0.01;
    

    // CONNECTION BETWEEN TOP AND BOTTOM (SIDE POLYGONS)

    double x1, y1, x2, y2, x1p = 0, y1p = 0, x2p = 0, y2p = 0; // current and previous points
    double x0 = 0, y0 = 0, x0p = 0 , y0p = 0; // real coordinates (without zoom)

    double x0c, y0c, x1c, x2c, y1c, y2c;
    // step through all corner points and calculate polygons

    if (tallness > 0)
    for (int cornern=-1;cornern<ncorners;cornern++){
      int cornern_ = cornern; if (cornern == -1) cornern_ = ncorners-1; // -1-th corner as the last corner
      array<double>@ r = o[cornern_];
      double radius = r[OC_RADIUS], ex = r[OC_CORNER_X], ey = r[OC_CORNER_Y];
      double xc = r[OC_CIRC_CENTER_X], yc = r[OC_CIRC_CENTER_Y];

      double radiusBot = radius*bot_zoom; // slightly less for smoother round corners
      double radiusTop = radius*top_zoom;
      bool is_rounded = (radius > 0);

      
      if (!is_rounded) {
        // for straight corners
        x0c = ex; y0c = ey;
        x1c = ex*bot_zoom; y1c = ey*bot_zoom;
        x2c = ex*top_zoom; y2c = ey*top_zoom;
      } else {
        // for rounded corners
        x0c = xc; y0c = yc;
        x1c = xc*bot_zoom; y1c = yc*bot_zoom;
        x2c = xc*top_zoom; y2c = yc*top_zoom;
      }
      // add rotation if needed
      if (is_rotated) {
        double v0 = sqrt(x0c*x0c + y0c*y0c); // real coordinates
        double v1 = v0*bot_zoom; //sqrt(x1c*x1c + y1c*y1c); // bottom surface
        double v2 = v0*top_zoom; // sqrt(x2c*x2c + y2c*y2c); // top surface
        double a1 = atan2(x0c,y0c);
        double si = sin(a1-rotR), co = cos(a1-rotR);
        x0c = si*v0; y0c = co*v0;
        x1c = si*v1; y1c = co*v1;
        x2c = si*v2; y2c = co*v2;
      }
      // add center offset
      x1c += bot_x_offset; y1c += bot_y_offset;
      x2c += top_x_offset; y2c += top_y_offset;
      // if (n == 1) drawPoint(x1c, y1c, 2);
      //drawPoint(x2c, y2c, 3, "#0000DD");

      uint points_on_round_corner = uint(radiusTop / 1.0); // detalizaton of round corners (how many polygons)
      if (points_on_round_corner < 3) points_on_round_corner = 3;
      if (points_on_round_corner > 30) points_on_round_corner = 30;
      // if (angle_toR < angle_fromR) angle_toR += twopi;

      // calculations for rounded corners
      double angle_fromR, angle_toR, angleDeltaR, curAngleR;
      if (is_rounded) {
          // now we know center points, lets calc points on radius
          angle_fromR = (r[OC_ARCANGLE_FROM_RAD] -pi2)*-1; // convert from Arc style angles
          angle_toR = (r[OC_ARCANGLE_TO_RAD] -pi2)*-1;
          angleDeltaR = (angle_fromR - angle_toR) / double(points_on_round_corner-1);
          curAngleR = angle_fromR;
      } else {
        // for straight corners
        points_on_round_corner = 1;
        x0 = x0c; y0 = y0c; x1 = x1c; y1 = y1c; x2 = x2c; y2 = y2c; 
      }
      
      // parsing round corner point and creating polygons
      // for straight corners pass once like there's one point in the round corner
      for (uint pointn = 0; pointn<points_on_round_corner; pointn++){
        if (is_rounded) {
          double si = sin(curAngleR-rotR), co = cos(curAngleR-rotR);
          x0 = x0c + si*radius;
          y0 = y0c + co*radius;
          x1 = x1c + si*radiusBot;
          y1 = y1c + co*radiusBot;
          x2 = x2c + si*radiusTop;
          y2 = y2c + co*radiusTop;
        }

        // adding a new polygon
        if ((pointn > 0) or ((cornern>-1) and !is_perfect_circle))  {
          // calc polygon properties
          bool polygon_visible = true;

          // normal (90 deg) angle to surface
          double angleNormR = atan2(x1p-x1,-(y1p-y1))+pi2;
          if (angleNormR > pi) angleNormR -= twopi;

          // calculate angle between top and bottom points
          // compare no normal and thus we can decide if polygon is visible
          double adiff1 = angDiffR(angleNormR, atan2(x1-x2,-(y1-y2)));
          //double angle2 = atan2(x2p-x1p,y2p-y1p);
          //double adiff2 = angDiffR(angleNormR, angle2);
          if (abs(adiff1)>pi2) polygon_visible = false; 

          if (polygon_visible) {
            // calculate polygon area
            //double pw = sqrt(pow((x0-x0p),2) + pow((y0-y0p),2));
            //double ph = tallness;
            // double ps = pw*ph; // area of the polygon

            // center of the polygon
            double pxc = x0p + (x0-x0p)/2, pyc = y0p + (y0-y0p)/2, pzc = (tallness)/2;

            // calculate Light properties for polygon
            double lightness = calcPolyLightnessRad(gui_x+pxc, gui_y+pyc, elevation+pzc, angleNormR, 0, bodyL, light_dependency);
            // status("Light intensity: " + Lintensity);

            //if (polygons.length==8) {
              //hue = 150;
              // status(" totalLight: " + totalLight);

            //   // the luminance of a polygon is as less, as area of polygon is less than normal area

            //   // Polygon: " + pw + " x " + ph + " = " + ps + 
            //   // status(" PL_angle_diffR: " + (PL_angle_diffR*180/pi) + " PL_alt_diffR: " + (PL_alt_diffR*180/pi) +  " Lintensity: " + Lintensity + " Ldistance: " + Ldistance + " Langle: " + (LangleR*180/pi) + " angleNorm: " + (angleNormR*180/pi));
            //   // drawPoint(pxc, pyc, 2, "#0000DD");
            //}

            convertHSLtoRGB(bodyH, bodyS, lightness, cr,cg,cb);
            
            visibility = 1;
            if (!polygon_visible) visibility = 0;

            // add polygon
            polygons.insertLast({POLY_TYPE_NORMAL, visibility, cr, cg, cb, ca, x1,y1,x2,y2,x2p,y2p,x1p,y1p, angleNormR});
          }
          if (cornern == (ncorners-1)) break; // we already done last corner when cornern was 0
        }
        
        // remember current points as "previous"
        x0p = x0; y0p = y0; x1p = x1; y1p = y1; x2p = x2; y2p = y2;

        // for rounded corners go to next angle
        if (is_rounded) curAngleR -= angleDeltaR;
      }
        
    } // end of corners loop

    // add top surface polygon
    // calculate Light properties for top polygon
    double hue = bodyH, sat = bodyS, orig_lightess = bodyL;
    // double lightIntensity = calcPolygonLightIntensity(gui_x, gui_y, elevation+tallness, 0, pi2);

    double lightness = calcPolyLightnessRad(gui_x, gui_y, tallelev, 0, pi2, orig_lightess, light_dependency*top_light_dependency_ratio);
    convertHSLtoRGB(hue, sat, lightness, cr,cg,cb);

    visibility = 1;

    // add decorated bevel
    if (has_decor_bevel) {
      polygons.insertLast({POLY_TYPE_DECOR_BEVEL, visibility, cr, cg, cb, ca, top_x_offset, top_y_offset, top_zoom, rotation});
    } 

    // add top polygon
    double top_zoom_minus_bevel = top_zoom;
    if (has_decor_bevel) top_zoom_minus_bevel *= (1-decorBevelWidth);
    if (top_zoom_minus_bevel > 0) {
      polygons.insertLast({POLY_TYPE_OBJECT2D, visibility, cr, cg, cb, ca, top_x_offset, top_y_offset, top_zoom_minus_bevel, rotation});
    }
    

  } // end of calcPolygons

  // calculate lightness [(hs)L] param for polygon
  double calcPolyLightnessRad(double x, double y, double  z, double angleNormR, double altNormR, double orig_lightness, double light_dependency_ratio = 0.2) {
    double lightIntensity = calcPolygonLightIntensity(x, y, z, angleNormR, altNormR);
    lightIntensity += Light.ambient_intensity;
    double lig = orig_lightness*pow(lightIntensity, light_dependency_ratio);
    if (lig > 1) lig = 1;
    //if (lig < 0) lig = 0;
    return lig;
  }

  // calculate lightness [(hs)L] param for polygon with angles in degrees
  double calcPolyLightness(double x, double y, double  z, double angleNorm, double altNorm, double orig_lightness, double light_dependency_ratio = 0.2) {
    double angleNormR = angleNorm*pi/180;
    double altNormR = altNorm*pi/180;
    return calcPolyLightnessRad(x, y,  z, angleNormR, altNormR, orig_lightness, light_dependency_ratio);
  }

  // calculate iluminance of polygon
  double calcPolygonLightIntensity(double x, double  y, double  z, double angleNormR, double altNormR, int light_source_n = -1){

    double LangleR, LaltR, Ldistance, Lsize, Lintensity, Lx, Ly, Lz;
    
    // add ambient light intensity
    double totalIntensity = 0;

    int nlights = int(Light.sources.length);
    int ls_start = 0;

    if (light_source_n > -1) {
      ls_start = light_source_n;
      if (nlights > light_source_n) nlights = light_source_n;
    }

    // step though all light sources
    for(int ls_n = ls_start; ls_n < nlights; ls_n++) {
      array<double>@ ls = Light.sources[ls_n];
      if (ls[LS_ENABLED] == 0) continue; // skip if light source is off

      // calculate LangleR, LaltR, Ldistance, Lsize, Lintensity, Lx, Ly, Lz for current light source and polygon position
      calcLight(x, y, z, ls, LangleR, LaltR, Ldistance, Lsize, Lintensity, Lx, Ly, Lz);

      // calculate 3d vectors from polygon (normal) and from polygon to light
      double pvi=0, pvj=0, pvk=0, lvi, lvj, lvk, PxQ, PQ, plv_angleR;
      
      // calculate normal vector from polygon (in 3d space)
      pvk=de(sin(altNormR));
      pvi=de(sin(angleNormR))*(1-pvk);
      pvj=-de(cos(angleNormR))*(1-pvk);
      if (altNormR > pi2) {
        pvi = -pvi; pvj = -pvj;
      }
      
      // this is vector to light source (in 3d space)
      lvi = de(Lx); lvj = de(Ly); lvk = de(Lz); 

      // angle between vectors in 3d space
      PxQ = (pvi*lvi) + (pvj*lvj) + (pvk*lvk);
      PQ = sqrt(pvi*pvi + pvj*pvj + pvk*pvk) * sqrt(lvi*lvi + lvj*lvj + lvk*lvk);
      plv_angleR = acos(PxQ/PQ);

      // amount of light falling on polygon depends on this 3d angle
      double lum_factor = cos(plv_angleR);
      if (lum_factor < 0) lum_factor = 0;
      totalIntensity += lum_factor*Lintensity;

      // how polygon is facing the light (how far angle is from normal 90 deg)
      // double PL_angle_diffR = angDiffR(LangleR, angleNormR);
      // double PL_alt_diffR = angDiffR(altNormR, LaltR);

      // double phv, pwv; 
      // // how much visible "width" and "height" changes to this light angle
      // if (altNormR == pi2) {
      //   // for polygons looking straight up (flat top surface for example)
      //   pwv = 1;
      //   phv = cos(PL_alt_diffR);
      // } else if (altNormR == 0) {
      //   // for polygons looking on sides (vertical, side polygons)
      //   pwv = cos(PL_angle_diffR);
      //   phv = cos(PL_alt_diffR);
      // }

      // if (pwv < 0) pwv = 0;
      // if (phv < 0) phv = 0;
      //double psv = pwv*phv; // area of polygon "visible" from current light angle

      // double lum_factor = psv;
      // double lum_factor = 0;
      // if (altNormR == pi2) { // facing straigh up
      //   lum_factor = cos(PL_alt_diffR);
      // } else {
      //   // this work for sides
      //   lum_factor = cos(PL_alt_diffR)*cos(PL_angle_diffR);
      // }

      // lum_factor = cos(PL_angle_diffR);
      // lum_factor = cos(PL_alt_diffR);

      // status(" PL_alt_diffR: " + PL_alt_diffR*180/pi + " lum_factor: " + lum_factor);

      // string text = "Lintensity: " + Lintensity;
      // text += " Ldistance: " + Ldistance;
      // text += "LangleR: " + (LangleR*180/pi);
      // text += " angleNormR: " + (angleNormR*180/pi);
      // text += " PL_angle_diffR: " + (PL_angle_diffR*180/pi);
      // text += " pw : " + pw + " pwv: " + pwv;
      // text += "LaltR: " + (LaltR*180/pi);
      // text += " altNormR: " + (altNormR*180/pi);
      // text += " PL_alt_diffR: " + (PL_alt_diffR*180/pi);
      // text += " ph : " + ph + " phv: " + phv;
      //text += " totalLight : " + totalLight;
      //status(text);
    }

    

    return totalIntensity;
  }

  //////////////////////////////////////////
  // DRAW METHODS 
  //////////////////////////////////////////

  // draw calculated shadows
  void drawShadows(Kt::Graphics::Context@ ctx){

    if (!has_shadows) return;

    // parse through all pre-calculated shadows
    int shadows_len = int(shadowsPreCalc.length);
    if (shadows_len == 0) return;

    if (is_rotated) ctx.transform.Rotate(rotation); 


    // prepare vars
    bool draw_only_top_shadow;
    int max_layers;
    double layer_opacity, single_layer_opacity, zoom_top_delta, zoom_top, zoom_bottom_delta, zoom_bottom, s1ratio, s2ratio;
    double s1xc, s1yc, s2xc, s2yc;
    double r, g, b;
    double start_op, end_op, op_delta;
    double p1s1x, p1s1y, p1s2x, p1s2y, p2s1x, p2s1y, p2s2x, p2s2y;
    


    for (int sh_n = 0; sh_n < shadows_len; sh_n++) {
      
      // array<double>@ sh = shadows[sh_n];
      array<double>@ s = shadowsPreCalc[sh_n];

      // prepare vars
      r = s[SC_COLOR_R]; g = s[SC_COLOR_G]; b = s[SC_COLOR_B];
      
      // prepare shadow path
      ctx.path.Clear();

      // check if shadows fully overlap
      draw_only_top_shadow = (s[SC_FULLY_OVERLAP] != 0);

      
      // draw several layers of shadow to emulate blur
      max_layers = int(s[SC_MAXLAYERS]);

      // opacity for a single layer
      single_layer_opacity = 1.0*s[SC_OPACITY]/double(max_layers);

      // opacity for current layer
      start_op = s[SC_GFROMOP]; end_op = s[SC_GTOOP];
      layer_opacity = start_op;
      op_delta = (end_op - start_op)/double(max_layers);

      // set initial zoom
      zoom_top_delta = s[SC_BLURZOOMDELTA];
      zoom_top = 1.0 - (double(max_layers)*zoom_top_delta*0.5);
      zoom_bottom_delta = s[SC_BLURZOOMDELTA]*elev_percent;
      zoom_bottom = 1.0 - (double(max_layers)*zoom_bottom_delta*0.5);


      // get shadow coords in local vars for faster access
      s1ratio = s[SC_S1_RATIO]; s1xc = s[SC_S1XC]; s1yc = s[SC_S1YC];
      s2ratio = s[SC_S2_RATIO]; s2xc = s[SC_S2XC]; s2yc = s[SC_S2YC];
      

      // read connection points
      p1s1x = s[SC_P1S1X]; p1s1y = s[SC_P1S1Y];
      p1s2x = s[SC_P1S2X]; p1s2y = s[SC_P1S2Y];
      p2s1x = s[SC_P2S1X]; p2s1y = s[SC_P2S1Y];
      p2s2x = s[SC_P2S2X]; p2s2y = s[SC_P2S2Y];

      // ctx.path.MoveTo(s1xc+30, s1yc+30);
      // ctx.source.SetRGBA(0,0,0,1);
      // ctx.WriteText("layers:" + max_layers );

      // status("elev_percent: " + elev_percent);
      for(int layer=0; layer < max_layers; layer++) {

        ctx.path.Clear();

        // if shadows don't fully overlap
        if (!draw_only_top_shadow) {
          // draw shadow 1 (bottom of object)
          drawObjectPath(ctx, s1xc*zoom_bottom, s1yc*zoom_bottom, s1ratio*zoom_bottom);
          //ctx.path.Close();

          // draw the connection between shadows
          ctx.path.NewSubPath();
          ctx.path.MoveTo(p1s1x*zoom_bottom, p1s1y*zoom_bottom);
          ctx.path.LineTo(p1s2x*zoom_top, p1s2y*zoom_top);
          ctx.path.LineTo(p2s2x*zoom_top, p2s2y*zoom_top);
          ctx.path.LineTo(p2s1x*zoom_bottom, p2s1y*zoom_bottom);
          //ctx.path.Close();
        } 

        // draw shadow 2 (top of object)
        ctx.path.NewSubPath();
        drawObjectPath(ctx, s2xc*zoom_top, s2yc*zoom_top, s2ratio*zoom_top);
        // ctx.path.Close();

        // set shadow source as simple
        ctx.source.SetRGBA(r, g, b, layer_opacity*single_layer_opacity);

        // set shadow fill as gradient
        // gradients eat significantly more CPU (if many layers especcially)
        // Kt::Graphics::GradientDrawPattern@ gradient;
        // @gradient = ctx.patterns.NewLinearGradient(s[SC_GFROMX]*zoom, s[SC_GFROMY]*zoom, s[SC_GTOX]*zoom, s[SC_GTOY]*zoom);
        // gradient.AddColorStopRGBA(0, s[SC_COLOR_R], s[SC_COLOR_G], s[SC_COLOR_B], s[SC_GFROMOP]*layer_opacity);
        // gradient.AddColorStopRGBA(1, s[SC_COLOR_R], s[SC_COLOR_G], s[SC_COLOR_B], s[SC_GTOOP]*layer_opacity);
        // gradient.SelectAsSource(); 

        // fill shadow
        //ctx.StrokePath();
        ctx.FillPath();

        // drawPoint(p1s1x,p1s1y, 3, "#DD0000");
        // drawPoint(p2s1x,p2s1y, 3, "#00DD00");

        // increase zoom for next layer
        zoom_top += zoom_top_delta;
        zoom_bottom += zoom_bottom_delta;
        // change next layer opacity
        layer_opacity += op_delta;

      } // end of layers cycle
      

    } // end of draw shadow cycle


    if (is_rotated) ctx.transform.Rotate(-rotation); 
  } // end of DrawShadows

  // main function to "draw" current object
  // (add object path, without actually filling or stroking it)
  void drawObjectPath(Kt::Graphics::Context@ ctx, double offset_x = 0, double offset_y = 0, double zoom = 1, double rotate = 0, bool do_transform_back = true) {

    // check if we have to modify (scale, zoom, rotate) object
    bool do_offset = false, do_zoom = false, do_rotate = false;

    // set tranformation flags
    if ((offset_x != 0) or (offset_y != 0)) do_offset = true;
    if (zoom != 1) do_zoom = true;
    if (rotate != 0) do_rotate = true;

    // perform transformations
    if (do_offset) ctx.transform.Translate(offset_x,offset_y);
    if (do_zoom) ctx.transform.Scale(zoom,zoom);
    if (do_rotate) ctx.transform.Rotate(rotate); 

    // not "draw" the actual object
    switch (ncorners) {

      case 4: // normal object with 4 angles
        if (round == 0) {
          // straight rectangle
          // ctx.path.Rectangle(-w/2, -h/2, w, h);
          ctx.path.MoveTo(o1x, o1y);
          ctx.path.LineTo(o2x, o2y);
          ctx.path.LineTo(o3x, o3y);
          ctx.path.LineTo(o0x, o0y);

        } else if ((round == 1) and (w == h)) {
          // straight circle
          double ro = o[0][OC_RADIUS];
          ctx.path.Arc(0, 0, ro, -90, -90.001);
        } else {
          double ro = o[0][OC_RADIUS];
          // rect with rounded corners
          ctx.path.MoveTo(oc0x, o0y);
          ctx.path.Arc(oc1x, oc1y, ro, -90, 0);
          ctx.path.Arc(oc2x, oc2y, ro, 0, 90);
          ctx.path.Arc(oc3x, oc3y, ro, 90, 180);
          ctx.path.Arc(oc0x, oc0y, ro, 180, 270);
        }
        break;

      default:
        {
          // if object has not 4 angles, draw it via corner points
          
          // start with first corner
          // string s;
          double r, r1 = o[0][OC_RADIUS];
          double x0 = o[0][OC_CORNER_X], y0 = o[0][OC_CORNER_Y];
          if (r1 == 0) {
            ctx.path.MoveTo(x0, y0);
          } else {
            ctx.path.Arc(o[0][OC_CIRC_CENTER_X], o[0][OC_CIRC_CENTER_Y], r1, o[0][OC_ARCANGLE_TO], o[0][OC_ARCANGLE_TO]);
          }
          // and line through all corners
          for(int n=1;n<ncorners;n++) {
            r = o[n][OC_RADIUS];
            if (r == 0) {
              ctx.path.LineTo(o[n][OC_CORNER_X], o[n][OC_CORNER_Y]);
            } else {
              ctx.path.Arc(o[n][OC_CIRC_CENTER_X], o[n][OC_CIRC_CENTER_Y], r, o[n][OC_ARCANGLE_FROM], o[n][OC_ARCANGLE_TO]);
            }
            
            //s+=" " + x+","+y+" " ;
          }
          // finish where we've started
          if (r1 == 0) {
            ctx.path.LineTo(x0, y0);
          } else {
            ctx.path.Arc(o[0][OC_CIRC_CENTER_X], o[0][OC_CIRC_CENTER_Y], r1, o[0][OC_ARCANGLE_FROM], o[0][OC_ARCANGLE_TO]);
          }
          //status_custom_text = s;
        }
      }


    // if we need to transform context back
    if (do_transform_back) {
      if (do_rotate) ctx.transform.Rotate(-rotate); 
      if (do_zoom) ctx.transform.Scale(1/zoom,1/zoom);
      if (do_offset) ctx.transform.Translate(-offset_x,-offset_y);
    }
  } // end og drawObjectPath



  void drawObjectDebug(Kt::Graphics::Context@ ctx){
    // and line through all corners
    for(int n=0;n<ncorners;n++) {
      // double r, g, b, h, s, l;
      // h = 0+(n*360/ncorners); s = 1; l = 0.5;
      // convertHSLtoRGB(h, s, l, r, g, b); 
      // drawPoint(o[n][OC_CORNER_X], o[n][OC_CORNER_Y], 3, "", 1, r, g, b);
      // ctx.path.MoveTo(o[n][OC_CORNER_X] + 3, o[n][OC_CORNER_Y] - 3); 
      // ctx.WriteText(""+o[n][OC_CORNER_ANGLE_RAD]*180/pi);
    }
    

    // array<double>@ ls = Light.sources[0];


    // status("PV: " + pvi + "," + pvj + ","+ pvk + " | LV:" + lvi + ", " + lvj + ", " + lvk + " | plv_angle: " + plv_angleR*180/pi);

    // double v = 50;
    // drawPoint(pvi*v, pvj*v, 3);
    // v = 80;
    // drawPoint(lvi*v, lvj*v, 5, "#FFAAFF");
  }

  // view light source on a canvas (for debugging)
  void drawLightSources(Kt::Graphics::Context@ ctx){
  
    // get number of light sources
    array<array<double>>@ lights = Light.sources;
    int nlights = int(lights.length);

    double lx, ly, intensity, size, r, g, b;

    for(int ln=0;ln<nlights;ln++) {

      array<double>@ ls = lights[ln];

      if ((ls[LS_ENABLED] > 0) and (ls[LS_PREVIEW] > 0)) {

        ctx.path.Clear();
        // coords relative to current canvas pos
        // if light is at pos x=0,y=0 it should be at center of gui
        lx = ls[LS_X] - gui_canvas_x;
        ly = ls[LS_Y] - gui_canvas_y; 
        // status("lx: " + ls[LS_X] + " ly: " + ls[LS_Y]);
        r = ls[LS_COLOR_R]; g = ls[LS_COLOR_G]; b = ls[LS_COLOR_B];
        intensity = ls[LS_INTENSITY];
        size = ls[LS_SIZE];
        // draw spreading light circles
        double total_circles = intensity*5;
        for(int n=1;n<total_circles;n++) {
          double cintensity = (pow(n, 1.5) * intensity)/(total_circles*total_circles)/5;
          double csize = size/2 + (total_circles*3 - 3*n)/1.5;
          ctx.path.Arc(lx, ly, csize, 0.001, 0);
          ctx.source.SetRGBA (r,g,b, cintensity);
          ctx.FillPath();
          ctx.path.Clear();
        }

        // draw main circle
        ctx.path.Arc(lx, ly, size/2, 0.001, 0);
        ctx.source.SetRGBA (r,g,b, intensity);
        ctx.FillPath();

        // draw stroke around light object
        ctx.source.SetRGBA (0,0,0, 0.3);
        ctx.settings.set_lineWidth(1);
        ctx.StrokePath();
      }
    }

  } // end of drawLightSources

  // draw the object (generic function)
  void DrawGeneric(Kt::Graphics::Context@ ctx) {
    // we do all drawing around center (0,0) (to be able to easily rotate object)
    // set anchor point and rotation
    ctx.transform.Translate(cw2+xc_offset,ch2+yc_offset);
    // if (is_rotated) ctx.transform.Rotate(rotation); 

    // draw shadows is any
    if (has_shadows) drawShadows(ctx);

    // draw 3d or 2d object
    if (enable_3d) {
      
      // in 3D mode draw the polygons
      for(uint pn=0;pn<polygons.length;pn++){
        array<double>@ poly = polygons[pn];
        int ptype = int(poly[POLY_TYPE]);

        // depending on type of polygon
        switch(ptype) {
          case POLY_TYPE_OBJECT2D:
          case POLY_TYPE_DECOR_BEVEL:
          {
            // draw object top or bottom surface
            bool visible = (poly[POLY_VISIBLE] > 0);
            if (visible) {
              // draw path
              ctx.path.Clear();
              drawObjectPath(ctx, poly[POLY_XOFFSET], poly[POLY_YOFFSET], poly[POLY_ZOOM], poly[POLY_ROTATION]);
              ctx.path.Close();
              // fill path
              ctx.source.SetRGBA(poly[POLY_COLOR_R], poly[POLY_COLOR_G], poly[POLY_COLOR_B], poly[POLY_COLOR_A]);
                ctx.FillPath();

              if (ptype == POLY_TYPE_DECOR_BEVEL) {
                drawDecorBevel(ctx);
              }
              
            }
            break;
          }

          case POLY_TYPE_NORMAL: {
            // draw "normal" (side) polygons
            double x1 = poly[POLY_X1], y1 = poly[POLY_Y1], x2 = poly[POLY_X2], y2 = poly[POLY_Y2];
            double x3 = poly[POLY_X3], y3 = poly[POLY_Y3], x4 = poly[POLY_X4], y4 = poly[POLY_Y4];
            bool visible = (poly[POLY_VISIBLE] > 0);

            if (visible) {
              ctx.path.Clear();
              ctx.path.MoveTo(x1, y1); // bottom point
              ctx.path.LineTo(x2, y2); // top point
              ctx.path.LineTo(x3, y3); // top previous
              ctx.path.LineTo(x4, y4); // bottom previous
              //ctx.path.Close();

              // stroke lines
              // ctx.source.SetRGBA(borderR,borderG,borderB, borderA);
              // ctx.StrokePath();
              // get polygon color
              // convertHSLtoRGB(bodyH+float(n+1)*40, bodyS, bodyL, cr,cg,cb);
              // convertHSLtoRGB(bodyH, bodyS, bodyL*0.6*(float(n+1)/float(polygons.length)), r,g,b);
              ctx.source.SetRGBA(poly[POLY_COLOR_R], poly[POLY_COLOR_G], poly[POLY_COLOR_B], poly[POLY_COLOR_A]);
              ctx.FillPath();

              // double path_length = (abs(x2-x1)+abs(y2-y1)+abs(x4-x3)+abs(y4-y3));
              // bool dont_stroke = (path_length < 10);
              ctx.settings.set_lineWidth(0.5);
              ctx.StrokePath();
            }

            if (pn == 0) {

              //ctx.path.MoveTo(x1+10, y1+100); 
              //ctx.source.SetRGBA(0,0,0,1);
              //ctx.WriteText(""+poly[POLY_ANGLENORMALRAD]*180/pi);

              //drawPoint(x1, y1, 3);
              // drawPoint(x2, y2, 1, "#DD0000");
              // drawPoint(x3, y3, 1, "#0000DD");
              //drawPoint(x4, y4, 3, "#FFAA00");
            }
            break;
          }
          
        }
      }
      // status("Polygons: " + polygons.length);

      // if shading/highlighting is needed
      if (has_shading or has_highlight) drawShading(ctx);

      // stroke object body top border
      ctx.settings.set_lineWidth(borderWidth);
      ctx.source.SetRGBA (borderR,borderG,borderB, borderA);
      ctx.StrokePath();



    } else {
      // draw 2D simple flat object
      if (is_rotated) ctx.transform.Rotate(rotation); 
      ctx.path.Clear();
      drawObjectPath(ctx);
      ctx.path.Close();
      // fill object body
      if (bodyA > 0) {
        ctx.source.SetRGBA (bodyR,bodyG,bodyB, bodyA);
        ctx.FillPath();
      }
      // stroke object border
      ctx.source.SetRGBA (borderR,borderG,borderB, borderA);
      ctx.settings.set_lineWidth(borderWidth);
      ctx.StrokePath();

      // if shading/highlighting is needed
      if (has_decor_bevel) drawDecorBevel(ctx);
      if (has_shading or has_highlight) drawShading(ctx);

      if (is_rotated) ctx.transform.Rotate(-rotation); 
    }

    ctx.source.SetRGBA(0,0,1,1);
    ctx.path.MoveTo(0,0);
    // ctx.WriteText("" + gui_canvas_x + "," + gui_canvas_y);
    if (debug_mode == 3) drawObjectDebug(ctx);

    // restore anchor point and rotation
    // if (is_rotated) ctx.transform.Rotate(-rotation); 
    ctx.transform.Translate(-(cw2+xc_offset),-(ch2+yc_offset));
  }
  
  
  void drawDecorBevel(Kt::Graphics::Context@ ctx){

    // draw highlight (light)
    auto bodyHLGrad = ctx.patterns.NewLinearGradient(ShadeX,ShadeY,-ShadeX,-ShadeY);
    
    // status("CL_r: " + CL_r + " " + CL_g + " " + CL_b);
    bodyHLGrad.AddColorStopRGBA(decorBevelHLPosA, CL_r, CL_g, CL_b, decorBevelHLPosAOp);
    bodyHLGrad.AddColorStopRGBA(decorBevelHLPosB, CL_r, CL_g, CL_b, decorBevelHLPosBOp);
    bodyHLGrad.SelectAsSource(); 
    ctx.FillPath();

    // draw bevel shade
    ctx.settings.set_blendMode(Kt::Graphics::kDrawOpXor);
    auto bodyShadeGrad = ctx.patterns.NewLinearGradient(ShadeX,ShadeY,-ShadeX,-ShadeY);
    bodyShadeGrad.AddColorStopRGBA(decorBevelShPosA, CL_r, CL_g, CL_b, decorBevelShPosAOp);
    bodyShadeGrad.AddColorStopRGBA(decorBevelShPosB, CL_r, CL_g, CL_b, decorBevelShPosBOp);
    bodyShadeGrad.SelectAsSource(); 
    ctx.FillPath();
    ctx.settings.set_blendMode(Kt::Graphics::kDrawOpOver);
  }

  void drawShading(Kt::Graphics::Context@ ctx) {
    double sx = ShadeX, sy = ShadeY;

    if ((has_shading) or (has_highlight)) {
      ctx.path.Clear();
      drawObjectPath(ctx, top_x_offset, top_y_offset, top_zoom, rotation);
      ctx.path.Close();
    }
    
    // we can swap shadine for "concave" effect
    if (swap_shade_edges) {
      sx = -ShadeX; 
      sy = -ShadeY;
    }

    // draw shade (black)
    if (has_shading) {
      ctx.settings.set_blendMode(Kt::Graphics::kDrawOpXor);
      auto bodyShadeGrad = ctx.patterns.NewLinearGradient(sx,sy,-sx,-sy);
      bodyShadeGrad.AddColorStopRGBA(0, CL_r, CL_g, CL_b, ShadeOpStart);
      bodyShadeGrad.AddColorStopRGBA(1, CL_r, CL_g, CL_b, ShadeOpEnd);
      bodyShadeGrad.SelectAsSource(); 
      ctx.FillPath();
      ctx.settings.set_blendMode(Kt::Graphics::kDrawOpOver);
    }

    // draw highlight (light)
    if (has_highlight) {
      ctx.settings.set_blendMode(Kt::Graphics::kDrawOpOver);
      auto bodyHLGrad = ctx.patterns.NewLinearGradient(sx,sy,-sx,-sy);
      bodyHLGrad.AddColorStopRGBA(0, CL_r, CL_g, CL_b, HLOpStart);
      bodyHLGrad.AddColorStopRGBA(1, CL_r, CL_g, CL_b, HLOpEnd);
      bodyHLGrad.SelectAsSource(); 
      ctx.FillPath();
    }

    // radial shade
    // } else {
    //   auto bodyHLRadGrad = ctx.patterns.NewRadialGradient(0, 0, 2, 0, 0, maxradius*2);
    //   bodyHLRadGrad.AddColorStopRGBA(0, 1, 1, 1, HLOpStart);
    //   bodyHLRadGrad.AddColorStopRGBA(1, 1, 1, 1, HLOpEnd);
    //   bodyHLRadGrad.SelectAsSource(); 
    //   ctx.FillPath();
    // }

    // drawLine(ShadeX,ShadeY,-ShadeX,-ShadeY);
    //drawPoint(maxradius * sin(CL_AngleR), maxradius * cos(CL_AngleR));

    //ctx.source.SetRGBA(1,0.5,0.5,1);
    //ctx.WriteText(round(HLOpStart,3) + " | " + round(HLOpEnd,3));
  }



} // end of CanvasObject Class


  // some extras here for dev purposes
  void drawPoint(double xc, double yc, double radius = 3, string color = "#00AA00", double a = 1, double r = -1, double g = -1, double b = -1) {
    auto ctx = Kt::Graphics::GetCurrentContext();
    ctx.SaveState();
    if (r == -1) {
      r=double(parseInt(color.substr(1,2),16))/255.0;
      g=double(parseInt(color.substr(3,2),16))/255.0;
      b=double(parseInt(color.substr(5,2),16))/255.0;
    }
    ctx.source.SetRGBA (r,g,b, a);
    ctx.path.Clear();
    ctx.path.Arc(xc, yc, radius, 0.001, 0);
    ctx.FillPath();
    ctx.RestoreState();
    ctx.path.Clear();
  }

  void drawLine(double x1, double y1, double x2, double y2, double width = 2, string color = "#0000DD", double a =1 ){
    auto ctx = Kt::Graphics::GetCurrentContext();
    ctx.SaveState();
    double r=double(parseInt(color.substr(1,2),16))/255.0;
    double g=double(parseInt(color.substr(3,2),16))/255.0;
    double b=double(parseInt(color.substr(5,2),16))/255.0;
    ctx.source.SetRGBA (r,g,b, a);
    ctx.path.Clear();
    ctx.path.MoveTo(x1, y1);
    ctx.path.LineTo(x2, y2);
    ctx.settings.set_lineWidth(width);
    ctx.StrokePath();
    ctx.RestoreState();
    ctx.path.Clear();
  }

}

