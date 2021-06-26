namespace LM{


  enum joystickParams{
    RJ_GUI_X,
    RJ_GUI_Y,
    RJ_PAD_TOP_PX,
    RJ_PAD_BOTTOM_PX,
    RJ_PAD_LEFT_PX,
    RJ_PAD_RIGHT_PX,
    RJ_THUMB_COLOR, 
    RJ_THUMB_OPACITY,
    RJ_THUMB_TYPE,
    RJ_THUMB_SIZE, 
    RJ_THUMB_SIZE_B_RATIO,
    RJ_THUMB_TALLNESS_RATIO,
    RJ_THUMB_ELEVATION_RATIO,
    RJ_THUMB_LEG_SIZE_RATIO,
    RJ_THUMB_BORDER_COLOR,
    RJ_THUMB_BORDER_OPACITY,
    RJ_THUMB_BORDER_WIDTH,
    RJ_THUMB_ROUND,
    RJ_THUMB_ENABLE_3D,
    RJ_SHADOW_TYPE,
    RJ_SHADOW_COLOR,
    RJ_SHADOW_INTENSITY,
    RJ_SHADOW_LIGHT_SOURCE,
    RJ_SHADOW_MAX_BLUR_LAYERS
  }

  // a joystick object
  class Joystick : CanvasWidget {

    double xnval = 0.5, ynval = 0.5, znval = 0.5; // normalized value of params (from 0 to 1)

    // thumb params
    joystickThumb thumb; // thumb object
    int thumbType = 0;
    float thumb_x = -1, thumb_y = -1; // center-pos recalculated for nval params
    float thumb_x_prev = -1, thumb_y_prev = -1; // to remember

    // Joystick constructor
    Joystick(float cw, float ch, string params){
      // remember canvas size
      this.cw = cw;
      this.ch = ch;
      cw2=cw/2.0;
      ch2=ch/2.0;
      
      // parse params
      array<string> ar = params.split(";");

      // padding
      padl = f(ar[RJ_PAD_LEFT_PX]);
      padr = f(ar[RJ_PAD_RIGHT_PX]);
      padt = f(ar[RJ_PAD_TOP_PX]);
      padb = f(ar[RJ_PAD_BOTTOM_PX]);
      // full canvas width
      cfw = cw + padl + padr;
      cfh = ch + padt + padb;

      // position of canvas top left about center of gui
      if (ar[RJ_GUI_X] == "") gui_canvas_x = -cfw/2; else gui_canvas_x = f(ar[RJ_GUI_X]);
      if (ar[RJ_GUI_Y] == "") gui_canvas_y = -cfh/2; else gui_canvas_y = f(ar[RJ_GUI_Y]);

      // thumb type
      thumbType = parseInt(ar[RJ_THUMB_TYPE]);
      // add thumb object if needed
      if (thumbType > 0) {
        // thumb params
        float thumb_size_b, thumb_border_w;
        float thumb_size = f(ar[RJ_THUMB_SIZE]);
        // border width depends on main size
        float thumb_border_w_ratio = f(ar[RJ_THUMB_BORDER_WIDTH]);
        thumb_border_w = thumb_border_w_ratio*thumb_size/10;
        // size b depends on main size
        thumb_size_b = thumb_size*f(ar[RJ_THUMB_SIZE_B_RATIO]);
        // tallness and elevation depends on main size
        float thumb_tallness = thumb_size*f(ar[RJ_THUMB_TALLNESS_RATIO]);
        float thumb_elevation = thumb_size*f(ar[RJ_THUMB_ELEVATION_RATIO]);

        bool thumb_enable_3d = (ar[RJ_THUMB_ENABLE_3D] == "1");
        // thumb_enable_3d = false;

        float thumb_rotation = 0;
        // initialize thumb
        thumb.init(cw, ch, gui_canvas_x, gui_canvas_y, thumb_enable_3d);
        thumb.setShape(4, thumb_size, thumb_size_b, thumb_tallness, thumb_elevation, f(ar[RJ_THUMB_ROUND]), thumb_rotation);
        thumb.setBodyColor(ar[RJ_THUMB_COLOR], f(ar[RJ_THUMB_OPACITY]));
        thumb.setBorder(ar[RJ_THUMB_BORDER_COLOR], f(ar[RJ_THUMB_BORDER_OPACITY]), thumb_border_w);
        //thumb.setMarker(thumb_size*f(ar[RJ_THUMB_MARKER_LENGTH_RATIO]), (thumb_size_b/10)*f(ar[RJ_THUMB_MARKER_THICKNESS_RATIO]), f(ar[RJ_THUMB_MARKER_ROUND]), ar[RJ_THUMB_MARKER_COLOR], f(ar[RJ_THUMB_MARKER_OPACITY]), "#000000", 0);
        //thumb.setDecor(parseInt(ar[RJ_THUMB_DECOR_TYPE]), f(ar[RJ_THUMB_DECOR_OPACITY]));
        //thumb.setShadeAndHighlight(f(ar[RJ_THUMB_SHADE_INTENSITY]));

        // add shadow (shadow type, color, opacity, supported light source (-1 for all sources), max_blur_layers)
        if ((ar[RJ_SHADOW_INTENSITY] != "0")) {
          string lss = ar[RJ_SHADOW_LIGHT_SOURCE];
          if (lss == "") lss = "-1"; // all light sources if not selected
          array<string> arlss = lss.split(",");
          for(uint sn=0;sn<arlss.length;sn++) {
            thumb.addShadow(parseInt(ar[RJ_SHADOW_TYPE]), ar[RJ_SHADOW_COLOR], f(ar[RJ_SHADOW_INTENSITY]), parseInt(arlss[sn]), parseInt(ar[RJ_SHADOW_MAX_BLUR_LAYERS]));
          }
        } 
          
      }


      // set nval to unknown
      updateNVals(UNKNOWN_VALUE, UNKNOWN_VALUE);
    }

    void updateNVals(double xnvalue, double ynvalue = 0, double znvalue = 0){
      this.xnval = xnvalue;
      this.ynval = ynvalue;
      this.znval = znvalue;
      
      if (xnval == UNKNOWN_VALUE) return; // xnval is not set ()

      thumb_x = xnval*cw + padl;
      thumb_y = ch-ynval*ch + padt;

      // set thumb center
      if ((thumb_y != thumb_y_prev) || (thumb_x != thumb_x_prev)) {
        thumb.setCenterPos(thumb_x, thumb_y);
      }

      thumb_x_prev = thumb_x; thumb_y_prev = thumb_y;

    }

    // on render settings change
    void updateRenderSettings() {
      thumb.updateRenderSettings();
    }

    // theoretically runs only once 
    void setInitValues(double initxnval = 0.5, double initynval = 0.5, double initznval = 0.5) {
      // if nval is not yet set (first render), set values
      if (xnval == UNKNOWN_VALUE) {
        updateNVals(initxnval, initynval, initznval);
      }
    }

    // main function to re-draw
    void Draw(Kt::Graphics::Context@ ctx){

      // draw the thumb
      thumb.Draw(ctx);

      // draw debug
      if (debug_mode != 0) DrawDebug(ctx);

    } // end of Joystick Draw
  } // end of Joystick class

  // joystickThumb extending CanvasObject
  class joystickThumb : CanvasObject {

    void Draw(Kt::Graphics::Context@ ctx) {
      DrawGeneric(ctx);
    }
  } 

} // end of LM namespace



