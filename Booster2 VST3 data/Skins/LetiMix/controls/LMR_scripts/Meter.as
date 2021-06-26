namespace LM{

  enum meterParams{
    RM_GUI_X,
    RM_GUI_Y,
    RM_GUI_XY_CENTERED,
    RM_PAD_TOP_PX,
    RM_PAD_BOTTOM_PX,
    RM_PAD_LEFT_PX,
    RM_PAD_RIGHT_PX
   
  }


  // a meter object
  class Meter : CanvasWidget {
    
    //double nval = 0.5; // normalized value of param (from 0 to 1)
    //double hold_nval = 0.5; // if hold is used

    // pre-calculated params
    int orientation = 0;

    // meter render params
    array<double> mp(80);

    // meter constructor
    Meter(float cw, float ch, int orient, string params, string render){
      // remember canvas size
      this.cw = cw;
      this.ch = ch;
      cw2=cw/2.0;
      ch2=ch/2.0;
      
      // parse params
      array<string> ar = params.split(";");

      // padding
      padl = f(ar[RM_PAD_LEFT_PX]);
      padr = f(ar[RM_PAD_RIGHT_PX]);
      padt = f(ar[RM_PAD_TOP_PX]);
      padb = f(ar[RM_PAD_BOTTOM_PX]);
      // full canvas width
      cfw = cw + padl + padr;
      cfh = ch + padt + padb;

      // position of meter's canvas top left about center of gui
      if (ar[RM_GUI_X] == "") gui_canvas_x = -cfw/2; else gui_canvas_x = f(ar[RM_GUI_X]);
      if (ar[RM_GUI_Y] == "") gui_canvas_y = -cfh/2; else gui_canvas_y = f(ar[RM_GUI_Y]);
      // get canvas left top position if it is given as center position
      if (ar[RK_GUI_XY_CENTERED] == "1") {
        gui_canvas_x -= cw2;
        gui_canvas_y -= ch2;
      }

      this.orientation = orient;

      // set meter render params
      string render_string = "$METERS_RENDER$";
      if (!render.isEmpty()) render_string = render;
      meters_prepareParams(render_string, mp);
      mp[RM_MIN_VAL] = 0;
      mp[RM_MAX_VAL] = 1;
    }


    // on light change
    void updateRenderSettings() {
      
    }

    // main function to re-draw the meter
    void Draw(Kt::Graphics::Context@ ctx, double nval, double hold_nval = UNKNOWN_VALUE){

      ctx.transform.Translate(padl,padt);
      renderMeter(ctx, cw, ch, orientation, mp, nval, hold_nval);

    } // end of Draw of Meter
  } // end of Meter class

}




