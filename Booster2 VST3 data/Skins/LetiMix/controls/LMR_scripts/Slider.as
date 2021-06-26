namespace LM{

  enum sliderParams{
    RS_GUI_X,
    RS_GUI_Y,
    RS_GUI_XY_CENTERED,
    RS_PAD_TOP_PX,
    RS_PAD_BOTTOM_PX,
    RS_PAD_LEFT_PX,
    RS_PAD_RIGHT_PX,
    RS_HORIZONTAL,
    RS_REVERSE,
    RS_SELECTION_COLOR,
    RS_SELECTION_OPACITY,
    RS_SELECTION_BG_COLOR,
    RS_SELECTION_BG_OPACITY,
    RS_SELECTION_WIDTH,
    RS_SELECTION_WIDTH_PX,
    RS_SELECTION_W_OFFSET,
    RS_SELECTION_LENGTH,
    RS_SELECTION_L_OFFSET,
    RS_SELECTION_ROUND,
    RS_SELECTION_SYM,
    RS_SELECTION_BORDER_WIDTH_RATIO,
    RS_SELECTION_BORDER_COLOR,
    RS_SELECTION_BORDER_OPACITY,
    RS_MARGIN_BEFORE_PX,
    RS_MARGIN_AFTER_PX,
    RS_THUMB_COLOR, 
    RS_THUMB_OPACITY,
    RS_THUMB_TYPE,
    RS_THUMB_SIZE_PX, 
    RS_THUMB_SIZE_B_RATIO,
    RS_THUMB_TALLNESS_RATIO,
    RS_THUMB_ELEVATION_RATIO,
    RS_THUMB_LEG_SIZE_RATIO,
    RS_THUMB_BORDER_COLOR,
    RS_THUMB_BORDER_OPACITY,
    RS_THUMB_BORDER_WIDTH,
    RS_THUMB_ROUND,
    RS_THUMB_MARKER_COLOR,
    RS_THUMB_MARKER_OPACITY,
    RS_THUMB_MARKER_LENGTH_RATIO,
    RS_THUMB_MARKER_THICKNESS_RATIO,
    RS_THUMB_MARKER_ROUND,
    RS_THUMB_DECOR_TYPE,
    RS_THUMB_DECOR_OPACITY,
    RS_THUMB_SHADE_INTENSITY,
    RS_THUMB_ENABLE_3D,
    RS_SHADOW_TYPE,
    RS_SHADOW_COLOR,
    RS_SHADOW_INTENSITY,
    RS_SHADOW_LIGHT_SOURCE,
    RS_SHADOW_MAX_BLUR_LAYERS,
    RS_USE_AS_BUTTON,
    RS_EXTRA_STYLE,
    RS_DEBUG_MODE
  }

  enum sliderThumbTypes{
    RS_TT_NONE,
    RS_TT_NORMAL,
    RS_TT_ONLY_BORDER 
  }

  // a slider object
  class Slider : CanvasWidget {
    
    double nval = 0.5; // normalized value of param (from 0 to 1)
    array<string> ar;

    // pre-calculated params
    bool is_vertical, is_horizontal, is_reversed, is_symmetric, has_selection_border = false, is_button = false, has_extra_style;
    string extra_style = "";
    float margin_before, margin_after;

    // selection params
    float sel_round, sel_width_ratio, sel_width_px, sel_w_offset, sel_length_ratio, sel_l_offset, selBgRound;
    float selbgx = 0, selbgy = 0, selbgw, selbgh; // normal selection size (useful)
    float selbgw_delta = 0, selbgh_delta = 0; // size with delta to compensate for thumb "leg"
    float dselbgx = 0, dselbgy = 0, dselbgw = 0, dselbgh = 0; // compensated for "leg" delta
    float selx, sely, selw, selh; // recalculated for nval params
    float selBorderWidth;

    float selBgR, selBgG, selBgB, selBgA; // colors
    float selR, selG, selB, selA; 
    float selBrdR, selBrdG, selBrdB, selBrdA;
    
    // thumb params
    SliderThumb thumb; // thumb object
    int thumbType = 0;
    float thumb_x = -1, thumb_y = -1; // center-pos recalculated for nval params
    float thumb_x_prev = -1, thumb_y_prev = -1; // to remember

    // slider constructor
    Slider(float cw, float ch, string params){
      // remember canvas size
      this.cw = cw;
      this.ch = ch;
      cw2=cw/2.0;
      ch2=ch/2.0;
      
      // parse params
      ar = params.split(";");

      // padding
      padl = f(ar[RS_PAD_LEFT_PX]);
      padr = f(ar[RS_PAD_RIGHT_PX]);
      padt = f(ar[RS_PAD_TOP_PX]);
      padb = f(ar[RS_PAD_BOTTOM_PX]);
      // full canvas width
      cfw = cw + padl + padr;
      cfh = ch + padt + padb;

      // position of slider's canvas top left about center of gui
      if (ar[RS_GUI_X] == "") gui_canvas_x = -cfw/2; else gui_canvas_x = f(ar[RS_GUI_X]);
      if (ar[RS_GUI_Y] == "") gui_canvas_y = -cfh/2; else gui_canvas_y = f(ar[RS_GUI_Y]);
      // get canvas left top position if it is given as center position
      if (ar[RK_GUI_XY_CENTERED] == "1") {
        gui_canvas_x -= cw2;
        gui_canvas_y -= ch2;
      }

      // general params
      is_vertical = (ar[RS_HORIZONTAL] == "0");
      is_horizontal = !is_vertical;
      is_reversed = (ar[RS_REVERSE] == "1");
      is_symmetric = (ar[RS_SELECTION_SYM] == "1");
      extra_style = ar[RS_EXTRA_STYLE];
      has_extra_style = (extra_style != "");
      

      // selection params
      margin_before = f(ar[RS_MARGIN_BEFORE_PX]);
      margin_after = f(ar[RS_MARGIN_AFTER_PX]);
      sel_round = f(ar[RS_SELECTION_ROUND]);
      sel_length_ratio = f(ar[RS_SELECTION_LENGTH]);
      sel_l_offset = f(ar[RS_SELECTION_L_OFFSET]);
      sel_w_offset = f(ar[RS_SELECTION_W_OFFSET]);
      sel_width_ratio = f(ar[RS_SELECTION_WIDTH]);
      sel_width_px = f(ar[RS_SELECTION_WIDTH_PX]);
      // if sel width is set in px, convert it to ratio
      if (sel_width_px > -0.5) {
        if (is_vertical) sel_width_ratio = sel_width_px/cw; 
        else sel_width_ratio = sel_width_px/ch;
      }

      // thumb params
      float thumb_size_b, thumb_border_w;
      float thumb_size = f(ar[RS_THUMB_SIZE_PX]);
      // border width depends on main size
      float thumb_border_w_ratio = f(ar[RS_THUMB_BORDER_WIDTH]);
      if (is_vertical) {
        //thumb_size = thumb_size_ratio*cw;
        thumb_border_w = thumb_border_w_ratio*thumb_size/10;
      } else {
        //thumb_size = thumb_size_ratio*ch;
        thumb_border_w = thumb_border_w_ratio*thumb_size/10;
      }
      // size b depends on main size
      thumb_size_b = thumb_size*f(ar[RS_THUMB_SIZE_B_RATIO]);
      // tallness and elevation depends on main size
      float thumb_tallness = thumb_size*f(ar[RS_THUMB_TALLNESS_RATIO]);
      float thumb_elevation = thumb_size*f(ar[RS_THUMB_ELEVATION_RATIO]);

      // thumb type
      thumbType = parseInt(ar[RS_THUMB_TYPE]);
      // add thumb object if needed
      if (thumbType > 0) {
        float thumb_rotation = 0;
        if (!is_vertical) thumb_rotation += 90;
        if (is_reversed) thumb_rotation -= 180;
        bool thumb_enable_3d = (ar[RS_THUMB_ENABLE_3D] == "1");

        // initialize thumb
        thumb.init(cw, ch, gui_canvas_x, gui_canvas_y, thumb_enable_3d);
        thumb.setShape(4, thumb_size, thumb_size_b, thumb_tallness, thumb_elevation, f(ar[RS_THUMB_ROUND]), thumb_rotation);
        thumb.setBodyColor(ar[RS_THUMB_COLOR], f(ar[RS_THUMB_OPACITY]));
        if (ar[RS_THUMB_BORDER_COLOR] == "") {
          ar[RS_THUMB_BORDER_COLOR] = ar[RS_SELECTION_COLOR];
          ar[RS_THUMB_BORDER_OPACITY] = ar[RS_SELECTION_OPACITY];
        }
        thumb.setBorder(ar[RS_THUMB_BORDER_COLOR], f(ar[RS_THUMB_BORDER_OPACITY]), thumb_border_w);
        thumb.setMarkerColors( ar[RS_THUMB_MARKER_COLOR], f(ar[RS_THUMB_MARKER_OPACITY]), "#000000", 0);
        thumb.setMarkerShapeForSlider(thumb_size*f(ar[RS_THUMB_MARKER_LENGTH_RATIO]), (thumb_size_b/10)*f(ar[RS_THUMB_MARKER_THICKNESS_RATIO]), f(ar[RS_THUMB_MARKER_ROUND]));
        thumb.setDecor(parseInt(ar[RS_THUMB_DECOR_TYPE]), f(ar[RS_THUMB_DECOR_OPACITY]));
        thumb.setShadeAndHighlight(f(ar[RS_THUMB_SHADE_INTENSITY]));

        // add shadow (shadow type, color, opacity, supported light source (-1 for all sources), max_blur_layers)
        if ((ar[RS_SHADOW_INTENSITY] != "0")) {
          string lss = ar[RS_SHADOW_LIGHT_SOURCE];
          if (lss == "") lss = "-1"; // all light sources if not selected
          array<string> arlss = lss.split(",");
          for(uint sn=0;sn<arlss.length;sn++) {
            thumb.addShadow(parseInt(ar[RS_SHADOW_TYPE]), ar[RS_SHADOW_COLOR], f(ar[RS_SHADOW_INTENSITY]), parseInt(arlss[sn]), parseInt(ar[RS_SHADOW_MAX_BLUR_LAYERS]));
          }
        } 
      }

      // selection area (whole under slider thumb)
      selbgw = cw; selbgh = ch;
      // calculate thumb "leg" size
      float thumb_leg_size = thumb_size_b*f(ar[RS_THUMB_LEG_SIZE_RATIO])*0.5;
      // calculate selection coordinates
      if (is_vertical) {
        // vertical slider
        selbgh = ch*sel_length_ratio - (margin_after + margin_before);
        selbgw = cw*sel_width_ratio;
        selbgx = (1-sel_width_ratio)*cw/2.0 + sel_w_offset*(cw-selbgw)/2.0 + padl;
        selbgy = (1-sel_length_ratio)*ch/2.0 + sel_l_offset*(ch-selbgh-(margin_after+margin_before))/2 + margin_after + padt;
        selbgh_delta = thumb_leg_size;
        dselbgx = selbgx; dselbgy = selbgy - selbgh_delta; dselbgw = selbgw; dselbgh = selbgh + selbgh_delta*2;
      } else {
        // horizontal slider
        selbgw = cw*sel_length_ratio - (margin_after + margin_before);
        selbgh = ch*sel_width_ratio;
        selbgy = (1-sel_width_ratio)*ch/2.0 + sel_w_offset*(ch-selbgh)/2.0 + padt;
        selbgx = (1-sel_length_ratio)*cw/2.0 + sel_l_offset*(cw-selbgw-(margin_after+margin_before))/2 + margin_before  + padl;
        selbgw_delta = thumb_leg_size;
        dselbgx = selbgx - selbgw_delta; dselbgy = selbgy; dselbgw = selbgw + selbgw_delta*2; dselbgh = selbgh;
      }
      selBgRound = sel_round;
      //if (selbgw > selbgh) selBgRound *= selbgh/2.0; else selBgRound *= selbgw/2.0;
      if (is_horizontal) {
        selBgRound *= selbgh/2.0; 
      } else {
        selBgRound *= selbgw/2.0;
      }
      //status("W" + selbgw + " H "+selbgh + " R " + selBgRound);
      
      if (is_horizontal) {
        selBorderWidth = selbgh*f(ar[RS_SELECTION_BORDER_WIDTH_RATIO]);
      } else {
        selBorderWidth = selbgw*f(ar[RS_SELECTION_BORDER_WIDTH_RATIO]);
      }

      // add padding to selection position
      // selbgx += padl; dselbgx += padl;
      //selbgy += padt; dselbgy += padt;


      // selection bg color
      selBgR=double(parseInt(ar[RS_SELECTION_BG_COLOR].substr(1,2),16))/255.0;
      selBgG=double(parseInt(ar[RS_SELECTION_BG_COLOR].substr(3,2),16))/255.0;
      selBgB=double(parseInt(ar[RS_SELECTION_BG_COLOR].substr(5,2),16))/255.0;
      selBgA=f(ar[RS_SELECTION_BG_OPACITY]);

      // selection color
      selR=double(parseInt(ar[RS_SELECTION_COLOR].substr(1,2),16))/255.0;
      selG=double(parseInt(ar[RS_SELECTION_COLOR].substr(3,2),16))/255.0;
      selB=double(parseInt(ar[RS_SELECTION_COLOR].substr(5,2),16))/255.0;
      selA=f(ar[RS_SELECTION_OPACITY]);

      // selection border color
      if (ar[RS_SELECTION_BORDER_COLOR] == "") {
        ar[RS_SELECTION_BORDER_COLOR] = ar[RS_SELECTION_COLOR];
        ar[RS_SELECTION_BORDER_OPACITY] = ar[RS_SELECTION_OPACITY];
      }
      selBrdR=double(parseInt(ar[RS_SELECTION_BORDER_COLOR].substr(1,2),16))/255.0;
      selBrdG=double(parseInt(ar[RS_SELECTION_BORDER_COLOR].substr(3,2),16))/255.0;
      selBrdB=double(parseInt(ar[RS_SELECTION_BORDER_COLOR].substr(5,2),16))/255.0;
      selBrdA=f(ar[RS_SELECTION_BORDER_OPACITY]);

      // calculate if we have selection border
      has_selection_border = ((selBorderWidth > 0) and (selBrdA > 0));
      
      // more flags
      is_button = (f(ar[RS_USE_AS_BUTTON]) > 0);

      debug_mode = parseInt(ar[RS_DEBUG_MODE]);
      // set nval to unknown
      updateNVal(UNKNOWN_VALUE);
    }

    // theoretically runs only once 
    void setInitValue(double initnval = 0.5) {
      // if nval is not yet set (first render), set values
      if (nval == UNKNOWN_VALUE) {
        updateNVal(initnval);
      }
    }

    // changes value of slider
    void updateNVal(double nvalue){
      this.nval = nvalue;
      
      if (nval == UNKNOWN_VALUE) return; // nval is not set ()
      if (nval<0) nval=0;
      if (nval>1) nval=1;

      // recalculate selection area
      selx = selbgx; sely = selbgy; selw = selbgw; selh = selbgh;
      if (is_vertical) {
        // vertical selection
        selx = selbgx; selw = selbgw;

        selh *= nval; 
        sely += (1-nval)*selbgh; 

        if (is_reversed) { 
          sely = selbgy;
          selh = (1-nval)*selbgh; 
        } 
      } else {
        // horizontal selection
        sely = selbgy; selh = selbgh;

        selw *= nval;
        if (is_reversed) {
          selx = selbgx + (nval)*selbgw;
          selw = selbgw - selx + selbgx; 
        } 
      }

      // calculate thumb center
      if (is_vertical) {
        thumb_x = selx + selw/2;
        if (!is_reversed) {
          thumb_y = sely;
        } else {
          thumb_y = sely+selh;
        }
      } else {
        thumb_y = sely + selh/2;
        if (is_reversed) {
          thumb_x = selx;
        } else {
          thumb_x = selx+selw;
        }
      }

      // set thumb center
      if ((thumb_y != thumb_y_prev) || (thumb_x != thumb_x_prev)) {
        thumb.setCenterPos(thumb_x, thumb_y);
        thumb.calcMarkerCoords();
        thumb.calcDecor();
      }

      thumb_x_prev = thumb_x; thumb_y_prev = thumb_y;

      // set thumb border color (for extra style)
      if (has_extra_style) {
        if (extra_style == "thumb_border_auto") {
          if (nval < 0.001) {
            thumb.borderR=double(parseInt(ar[RS_SELECTION_BG_COLOR].substr(1,2),16))/255.0;
            thumb.borderG=double(parseInt(ar[RS_SELECTION_BG_COLOR].substr(3,2),16))/255.0;
            thumb.borderB=double(parseInt(ar[RS_SELECTION_BG_COLOR].substr(5,2),16))/255.0;
            thumb.borderA=f(ar[RS_SELECTION_BG_OPACITY]);
          } else {
            thumb.borderR=double(parseInt(ar[RS_THUMB_BORDER_COLOR].substr(1,2),16))/255.0;
            thumb.borderG=double(parseInt(ar[RS_THUMB_BORDER_COLOR].substr(3,2),16))/255.0;
            thumb.borderB=double(parseInt(ar[RS_THUMB_BORDER_COLOR].substr(5,2),16))/255.0;
            thumb.borderA=f(ar[RS_THUMB_BORDER_OPACITY]);
          }
        }
      }
    }

    // on light change
    void updateRenderSettings() {
      thumb.updateRenderSettings();
      thumb.calcDecor();
    }

    void DrawSelectionPath(Kt::Graphics::Context@ ctx){
      ctx.path.Clear();
      if (selBgRound > 0) {
        // rounded selection
        float ro = selBgRound;
        ctx.path.ArcNegative(dselbgx+ro, dselbgy+ro, ro, -90, -180);
        ctx.path.ArcNegative(dselbgx+ro, dselbgy+dselbgh-ro, ro, 180, 90);
        ctx.path.ArcNegative(dselbgx+dselbgw-ro, dselbgy+dselbgh-ro, ro, 90, 0);
        ctx.path.ArcNegative(dselbgx+dselbgw-ro, dselbgy+ro, ro, 0, -90);
        ctx.path.Close();
      } else {
        // rectangular selection
        ctx.path.Rectangle(dselbgx, dselbgy, dselbgw, dselbgh);
      }
    }

    // main function to re-draw the slider
    void Draw(Kt::Graphics::Context@ ctx){

      // draw mask for selection bg
      DrawSelectionPath(ctx);

      // ctx.FillPath();
      
      // save state before clip with path
      ctx.SaveState();
      ctx.ClipWithPath(); // now we can draw only inside selection bg mask
      
      // PATCH: special "empty space" in mask for EMPTY-CENTERED thumbs
      if (thumbType == RS_TT_ONLY_BORDER) {
        ctx.path.Clear();
        float empty_space = thumb.h/2;
        if (is_vertical) {
          // for vertical sliders
          float h1 = selbgh*(1-nval)-empty_space;
          // if (is_reversed) h1 = selbgh*nval-empty_space; 
          h1 += selbgh_delta;
          if (h1>0) ctx.path.Rectangle(selbgx,selbgy-selbgh_delta, selbgw, h1);
          float h2 = selbgh-h1+empty_space*2;
          if (h2>0) ctx.path.Rectangle(selbgx,selbgy+h1+empty_space*2-selbgh_delta, selbgw, h2);
        } else {
          // for horizontal sliders
          float w1 = selbgw*nval-empty_space;
          //if (!is_reversed) w1 = selbgw*nval-empty_space; 
          w1 += selbgw_delta;
          if (w1>0) ctx.path.Rectangle(selbgx-selbgw_delta,selbgy, w1, selbgh);
          float w2 = selbgw-w1;
          if (w2>0) ctx.path.Rectangle(selbgx+w1+empty_space*2-selbgw_delta,selbgy,w2,selbgh);
        }
        ctx.ClipWithPath();
      }
      
      // draw selection bg
      ctx.path.Clear();
      ctx.source.SetRGBA (selBgR,selBgG,selBgB, selBgA);
      ctx.path.Rectangle(0, 0, cfw, cfh);
      ctx.FillPath();

      // draw selection
      ctx.path.Clear();
      ctx.source.SetRGBA (selR,selG,selB, selA);
      if (!is_symmetric) {
        // normal selection, add delta to compensate for thumb "leg"
        if (!is_reversed) { 
          ctx.path.Rectangle(selx-selbgw_delta, sely, selw+selbgw_delta, selh+selbgh_delta);
        } else {
          ctx.path.Rectangle(selx, sely-selbgh_delta, selw+selbgw_delta, selh+selbgh_delta);
        }
      } else {
        // symmetric selection
        if (is_vertical) {
          float selyc = selbgy+selbgh/2;
          float selhsym = selyc-selh-selbgy; if (is_reversed) selhsym = -selhsym;
          ctx.path.Rectangle(selx, selyc, selw, selhsym);
        } else {
          float selxc = selbgx+selbgw/2;
          float selwsym = selxc-selw-selbgx; if (!is_reversed) selwsym = -selwsym;
          ctx.path.Rectangle(selxc, sely, selwsym, selh);
        }
      }
      // draw selection
      ctx.FillPath();

      // undo clip masking
      ctx.RestoreState();

      // draw selection border
      if (has_selection_border) {
        DrawSelectionPath(ctx);
        //ctx.SaveState();
        //ctx.ClipWithPath();
        ctx.source.SetRGBA(selBrdR, selBrdG, selBrdB, selBrdA);
        ctx.settings.set_lineWidth(selBorderWidth);
        ctx.StrokePath();
        //ctx.RestoreState();
      }

      // draw the thumb
      thumb.Draw(ctx);

      // draw debug
      if (debug_mode != 0) DrawDebug(ctx);

    } // end of Draw of Slider
  } // end of Slider class


  // Slider Object extending CanvasObject
  class SliderThumb : CanvasObject {

    // marker for slider
    double markerThinkness = -1, markerLength = 0, markerRound = 0;
    double markerThinkness_orig, markerLength_orig;
    double mx0, my0, mx1, my1, mro; // marker edge points (calculated)
    
    // marker decoration
    double decorOpacity;
    int decorType;

    // gradients
    Kt::Graphics::GradientDrawPattern@ bodyDecorGrad, bodyShadeGrad, cur_grad, decorGradA, decorGradB, decorGradC;

    double grad_op_ratio_mod, grad_li_mid_point, grad_black_ratio, grad_white_ratio;
    double grad_max_white = 1, grad_max_black = -0.0;
    double li1, li2, li3, li4, lic; // light intensity of polygons

    // add marker information
    void setMarkerShapeForSlider(double markerLength, double markerThinkness, double markerRound) {
      
      this.markerThinkness_orig = markerThinkness;
      this.markerLength_orig = markerLength;

      this.markerThinkness = markerThinkness;
      this.markerLength = markerLength;
      this.markerRound = markerRound;
      calcMarkerCoords();
    }

    // add decor information
    void setDecor(int decor_type, double decor_opacity){
      decorType = decor_type;
      decorOpacity = decor_opacity;
    }

    void calcMarkerCoords(){
      if (!has_marker) return;

      // pre-calculate coordinates
      mx0 = (-markerLength/2);
      my0 = (-markerThinkness/2);
      mx1 = (markerLength/2);
      my1 = (markerThinkness/2);
      if (markerLength > markerThinkness) 
        mro = markerRound*markerThinkness/2.0; 
      else 
        mro = markerRound*markerLength/2.0;
    }
    
    // calculate params for gradients
    void calcDecor(){
      // lightness in center of gui
      lic = calcPolygonLightIntensity(0, 0, 0, 0, 90);

      // calculate lighness for "polygons" of fader
      li1 = calcLi(0, -h2, 0, 30); // edge top
      li2 = calcLi(0, -h2/2, 180, 30); // mid top
      li3 = calcLi(0, h2/2, 0, 30); // mid bot
      li4 = calcLi(0, h2, 180, 30); // bot
      // status("lic: " + lic);

    }

    ///////////////////////
    // FOR ADDING GRADIENTS

    // set current working gradient 
    void gradCurrent(Kt::Graphics::GradientDrawPattern@ gradient, double op_ratio_mod = 2, double li_mid_point = 0.23, double white_ratio = 1, double black_ratio = 1, double max_w = 1, double max_b = -1){
      @this.cur_grad = gradient;
      grad_op_ratio_mod = op_ratio_mod;
      grad_li_mid_point = li_mid_point;
      grad_white_ratio = white_ratio;
      grad_black_ratio = black_ratio;
      grad_max_white = max_w;
      grad_max_black = max_b;
    }

    // add gradient point and highlight (from -1 to 1)
    void grad(double p, double hl){
      double r = 1, g = 1, b = 1;
      if (hl > grad_max_white) hl = grad_max_white;
      if (hl < grad_max_black) hl = grad_max_black;
      if (hl<0) { 
        r=0; g=0; b=0;
      }
      cur_grad.AddColorStopRGBA(p, r, g, b, abs(hl)*grad_op_ratio_mod);
    }

    // add black or white starting point depending on lightness
    void grad0(double p, double lig){
      double hl = (lig-grad_li_mid_point)*0.0001;
      grad(p, hl);
    }

    // converting light intensity to grad lightness (-1..1)
    void grad2(double p, double lig, double mult = 1){
      double hl = lig-grad_li_mid_point;
      if (hl<0) {
        hl *= grad_black_ratio*mult;
      } else {
        hl *= grad_white_ratio*mult;
      }
      grad(p,hl);
    }

    // calculate lightness for given point on thumb and angle in degrees of virtual "polygon"
    // "polygon" xo, yo, angle, alt
    double calcLi(double xo, double yo, double angle, double alt, double topow = 1) {
      double z = elevation+tallness;
      // add rotation if needed
      if (is_rotated) {
        double v = sqrt(xo*xo + yo*yo); // vector
        double a1 = atan2(xo,yo);
        double si = sin(a1-rotR), co = cos(a1-rotR);
        xo = si*v; yo = co*v;
      }
      double x = gui_x + xo;
      double y = gui_y + yo;
      // drawPoint(x, y);      ctx.path.MoveTo(x,y);      ctx.WriteText("" + x + "," + y);
      double angleR = angle*pi/180;
      double altR = alt*pi/180;
      double li = calcPolygonLightIntensity(x, y, z, angleR+rotR, altR);
      if (topow != 1) li = pow(li, topow);
      return li;
    }

    void gradApply(Kt::Graphics::Context@ ctx) {
      cur_grad.SelectAsSource();
      ctx.FillPath();
    }

    // END OF GRADIENT RELATED STUFF
    //////////


    // draw the object
    void Draw(Kt::Graphics::Context@ ctx) {
      DrawGeneric(ctx);

      ctx.transform.Translate(cw2+xc_offset+top_x_offset,ch2+yc_offset+top_y_offset);
      ctx.transform.Rotate(rotation); 
      ctx.transform.Scale(top_zoom, top_zoom);

      // draw thumb decor
      bool clipping_on = false;
      if (decorType > 0) {
        ctx.SaveState();  ctx.ClipWithPath(); clipping_on = true;
        ctx.path.Clear();
        ctx.path.Rectangle(-cw2,-ch2,cw,ch);
        float aa; // quick alfa adjust
        switch(decorType) {

          case 1: // black slider 1
          // overall lighness
          @decorGradA = ctx.patterns.NewLinearGradient(o0x,o0y,o3x,o3y);
          gradCurrent(decorGradA, 1.8*decorOpacity, lic/1.4, 1, 0.4, 1, -.1);
          grad2(0.0, li1);
          grad2(0.13, li1, 0.8);
          grad2(0.20, li1);
          grad2(0.25, li1, 1.2);
          grad0(0.29, li1);
          grad2(0.29, lic, 0.1);
          grad2(0.71, lic, 0.1);
          grad0(1-0.28, li4);
          grad2(1-0.25, li4, 1.2);
          grad2(1-0.20, li4);
          grad2(1-0.13, li4, 0.8);     
          grad2(1.0, li4);
          gradApply(ctx);

          // additional shading
          @decorGradB = ctx.patterns.NewLinearGradient(o0x,o0y,o3x,o3y);
          gradCurrent(decorGradB, 0.65*decorOpacity); 
          grad(0.0, 0.1);
          grad(0.12, 0.09);
          grad(0.23, 0.15);
          grad(0.3, 0);
          grad(0.7, 0);
          grad(0.77, 0.15);
          grad(0.88, 0.09);
          grad(1, 0.1);
          gradApply(ctx);

          break;

          /* old style of type 1 (no lighting) */
          // aa=0.4*decorOpacity;
          // @bodyDecorGrad = ctx.patterns.NewLinearGradient(o0x,o0y,o3x,o3y);
          // bodyDecorGrad.AddColorStopRGBA(0.0, 1, 1, 1, 0.45*aa);
          // bodyDecorGrad.AddColorStopRGBA(0.13, 1, 1, 1, 0.4*aa);
          // bodyDecorGrad.AddColorStopRGBA(0.20, 1, 1, 1, 0.45*aa);
          // bodyDecorGrad.AddColorStopRGBA(0.25, 1, 1, 1, 0.6*aa);
          // bodyDecorGrad.AddColorStopRGBA(0.28, 1, 1, 1, 0*aa);

          // bodyDecorGrad.AddColorStopRGBA(0.75, 1, 1, 1, 0*aa);
          // bodyDecorGrad.AddColorStopRGBA(0.76, 1, 1, 1, 0.22*aa);
          // bodyDecorGrad.AddColorStopRGBA(0.89, 1, 1, 1, 0.33*aa);
          // bodyDecorGrad.AddColorStopRGBA(1, 1, 1, 1, 0.32*aa);
          // bodyDecorGrad.SelectAsSource(); 
          // ctx.FillPath();


          case 2:
          {
          // overall lighness
          @decorGradA = ctx.patterns.NewLinearGradient(o0x,o0y,o3x,o3y);
          gradCurrent(decorGradA, 5*decorOpacity, lic/1.4, 1, 0.4, 1, -.1);
          grad2(0.0, li1);
          grad2(0.12, li1);
          grad2(0.13, li2, 1.5);
          grad0(0.49, li2);
          grad0(0.51, li3);
          grad2(0.86, li3, 1.5);
          grad2(0.88, li4);
          grad2(1, li4);
          gradApply(ctx);

          // additional shading
          @decorGradB = ctx.patterns.NewLinearGradient(o0x,o0y,o3x,o3y);
          gradCurrent(decorGradB, 2*decorOpacity); 
          grad(0.0, -0.2);
          grad(0.1, -0.01);
          grad(0.14, 0.01);
          grad(0.5, 0.1);
          grad(0.86, 0.01);
          grad(0.9, -0.01);
          grad(1, -0.2);
          gradApply(ctx);

          // lines
          @decorGradC = ctx.patterns.NewLinearGradient(o0x,o0y,o3x,o3y);
          gradCurrent(decorGradC, 0.5*decorOpacity); 
          grad(0.185, -0.01); grad(0.195, -0.22); grad(0.205, -0.22); grad(0.215, -0.01);
          grad(0.275, -0.01); grad(0.28, -0.2); grad(0.30, -0.2); grad(0.32, -0.01);
          grad(0.365, -0.01); grad(0.385, -0.2); grad(0.405, -0.2); grad(0.41, -0.01);
          grad(0.46, -0.01); grad(0.47, -0.7); grad(0.53, -0.7); grad(0.54, -0.01);
          grad(0.59, -0.01); grad(0.595, -0.2); grad(0.615, -0.2); grad(0.635, -0.01);
          grad(0.685, -0.01); grad(0.705, -0.2); grad(0.725, -0.2); grad(0.73, -0.01);
          grad(0.785, -0.01); grad(0.795, -0.22); grad(0.805, -0.22); grad(0.815, -0.01);
          gradApply(ctx);
          
          }
          break;


          case 3:
          {
          // lines
          @decorGradC = ctx.patterns.NewLinearGradient(o0x,o0y,o3x,o3y);
          gradCurrent(decorGradC, 0.5*decorOpacity); 

          grad(0.23, 0.01); grad(0.24, 0.5); grad(0.26, 0.5); grad(0.27, 0.01);
          grad(0.47, 0.01); grad(0.48, 1); grad(0.52, 1); grad(0.53, 0.01);
          grad(0.73, 0.01); grad(0.74, 0.5); grad(0.76, 0.5); grad(0.77, 0.01);

          gradApply(ctx);
          
          }
          break;


          case 4: // metals
          {
          // @decorGradC = ctx.patterns.NewRadialGradient(o3x,o3y,20, o3x*0.1,o3y*0.1,1);
          // gradCurrent(decorGradC, 1*decorOpacity); 

          // grad(0.0, -0.01); grad(1, -0.7);


          // gradApply(ctx);
          
          }
          break;

          /*
          case 2: // old slider decor (light insensitive)

          aa=1*decorOpacity;
          @bodyDecorGrad = ctx.patterns.NewLinearGradient(o0x,o0y,o3x,o3y);
          bodyDecorGrad.AddColorStopRGBA(0.0, 0, 0, 0, 0.2*aa);
          bodyDecorGrad.AddColorStopRGBA(0.07, 1, 1, 1, 0.2*aa);
          bodyDecorGrad.AddColorStopRGBA(0.12, 1, 1, 1, 0.4*aa);
          bodyDecorGrad.AddColorStopRGBA(0.13, 0, 0, 0, 0.2*aa);

          bodyDecorGrad.AddColorStopRGBA(0.5, 1, 1, 1, 0.5*aa);

          bodyDecorGrad.AddColorStopRGBA(0.86, 1, 1, 1, 0.4*aa);
          bodyDecorGrad.AddColorStopRGBA(0.88, 0, 0, 0, 0.24*aa);
          bodyDecorGrad.AddColorStopRGBA(1, 0, 0, 0, 0.24*aa);
          bodyDecorGrad.SelectAsSource(); 
          ctx.FillPath();

          aa=0.3*decorOpacity;
          @bodyDecorGrad = ctx.patterns.NewLinearGradient(o0x,o0y,o3x,o3y);

          bodyDecorGrad.AddColorStopRGBA(0.185, 0, 0, 0, 0.0*aa);
          bodyDecorGrad.AddColorStopRGBA(0.195, 0, 0, 0, 0.22*aa);
          bodyDecorGrad.AddColorStopRGBA(0.205, 0, 0, 0, 0.22*aa);
          bodyDecorGrad.AddColorStopRGBA(0.215, 0, 0, 0, 0.0*aa);

          bodyDecorGrad.AddColorStopRGBA(0.275, 0, 0, 0, 0.0*aa);
          bodyDecorGrad.AddColorStopRGBA(0.28, 0, 0, 0, 0.2*aa);
          bodyDecorGrad.AddColorStopRGBA(0.30, 0, 0, 0, 0.2*aa);
          bodyDecorGrad.AddColorStopRGBA(0.32, 0, 0, 0, 0.0*aa);

          bodyDecorGrad.AddColorStopRGBA(0.365, 0, 0, 0, 0.0*aa);
          bodyDecorGrad.AddColorStopRGBA(0.385, 0, 0, 0, 0.2*aa);
          bodyDecorGrad.AddColorStopRGBA(0.405, 0, 0, 0, 0.2*aa);
          bodyDecorGrad.AddColorStopRGBA(0.41, 0, 0, 0, 0.0*aa);

          bodyDecorGrad.AddColorStopRGBA(0.46, 0, 0, 0, 0.0*aa);
          bodyDecorGrad.AddColorStopRGBA(0.47, 0, 0, 0, 0.7*aa);
          bodyDecorGrad.AddColorStopRGBA(0.53, 0, 0, 0, 0.7*aa);
          bodyDecorGrad.AddColorStopRGBA(0.54, 0, 0, 0, 0.0*aa);

          bodyDecorGrad.AddColorStopRGBA(0.59, 0, 0, 0, 0.0*aa);
          bodyDecorGrad.AddColorStopRGBA(0.595, 0, 0, 0, 0.2*aa);
          bodyDecorGrad.AddColorStopRGBA(0.615, 0, 0, 0, 0.2*aa);
          bodyDecorGrad.AddColorStopRGBA(0.635, 0, 0, 0, 0.0*aa);

          bodyDecorGrad.AddColorStopRGBA(0.685, 0, 0, 0, 0.0*aa);
          bodyDecorGrad.AddColorStopRGBA(0.705, 0, 0, 0, 0.2*aa);
          bodyDecorGrad.AddColorStopRGBA(0.725, 0, 0, 0, 0.2*aa);
          bodyDecorGrad.AddColorStopRGBA(0.73, 0, 0, 0, 0.0*aa);

          bodyDecorGrad.AddColorStopRGBA(0.785, 0, 0, 0, 0.0*aa);
          bodyDecorGrad.AddColorStopRGBA(0.795, 0, 0, 0, 0.22*aa);
          bodyDecorGrad.AddColorStopRGBA(0.805, 0, 0, 0, 0.22*aa);
          bodyDecorGrad.AddColorStopRGBA(0.815, 0, 0, 0, 0.0*aa);


          bodyDecorGrad.SelectAsSource(); 
          ctx.FillPath();

          break;
          */
        }
      }

      // draw marker
      if ((markerThinkness > 0) and (markerA > 0)) {
       ctx.path.Clear();
       ctx.source.SetRGBA (markerR,markerG,markerB, markerA);
       ctx.path.ArcNegative(mx0+mro, my0+mro, mro, -90, -180);
       ctx.path.ArcNegative(mx0+mro, my1-mro, mro, 180, 90);
       ctx.path.ArcNegative(mx1-mro, my1-mro, mro, 90, 0);
       ctx.path.ArcNegative(mx1-mro, my0+mro, mro, 0, -90);
       ctx.path.Close();
       ctx.FillPath();
      }

     // if shading is needed
     if (shadeIntensity > 0) {
       // if (!clipping_on) { ctx.SaveState(); ctx.ClipWithPath(); clipping_on = true; }
       // ctx.path.Clear();
       // ctx.path.Rectangle(-cw2,-ch2,cw,ch);
       // ctx.settings.set_blendMode(Kt::Graphics::kDrawOpXor);
     //   // calculate edge coordinates for shade
     //   float v = w2; if (h2 > w2) v = h2; // length of vector till longest edge
     //   float laRN = (LangleR+pi) - rotR; // light angle excluding rotation of thumb
     //   // calculate shading gradient edge points
     //   float ge1x, ge1y, ge4x, ge4y; // shading gradient edge points
     //   ge1x = sin(laRN)*v; ge1y = cos(laRN)*v; ge4x = sin(laRN+pi)*v; ge4y = cos(laRN+pi)*v;
     //   @bodyShadeGrad = ctx.patterns.NewLinearGradient(ge1x,ge1y,ge4x,ge4y);

     //   // reversed @bodyShadeGrad = ctx.patterns.NewLinearGradient(ge4x,ge4y,ge1x,ge1y);
     //   bodyShadeGrad.AddColorStopRGBA(0.1, 1, 1, 1, 0);
     //   bodyShadeGrad.AddColorStopRGBA(1, 1, 1, 1, shadeIntensity);
       // bodyShadeGrad.SelectAsSource(); 
       // ctx.FillPath();
     }

      // undo "clip with path" if on
      if (clipping_on) { ctx.RestoreState(); clipping_on = false; }

      ctx.transform.Scale(1/top_zoom, 1/top_zoom);
      ctx.transform.Rotate(-rotation); 
      ctx.transform.Translate(-(cw2+xc_offset+top_x_offset),-(ch2+yc_offset+top_y_offset));

      // for debugging
      if (true) drawLightSources(ctx);
    }


    

  } // end of class SliderThumb
}




