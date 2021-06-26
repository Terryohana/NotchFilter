namespace LM{

  enum knobParams{
    RK_GUI_X,
    RK_GUI_Y,
    RK_GUI_XY_CENTERED,
    RK_PAD_TOP_PX,
    RK_PAD_BOTTOM_PX,
    RK_PAD_LEFT_PX,
    RK_PAD_RIGHT_PX,

    RK_BODY_ENABLE_3D,
    RK_BODY_SIZE_RATIO,
    RK_BODY_TALLNESS_RATIO,
    RK_BODY_ELEVATION_RATIO,
    RK_BODY_COLOR,
    RK_BODY_OPACITY,
    RK_BODY_BORDER_COLOR,
    RK_BODY_BORDER_OPACITY,
    RK_BODY_BORDER_WIDTH,
    RK_BODY_CORNERS,
    RK_BODY_ROUND,

    RK_ANGLE_START, 
    RK_ANGLE_END,

    RK_MARKER_TYPE,
    RK_MARKER_COLOR, 
    RK_MARKER_OPACITY,
    RK_MARKER_SIZE_RATIO,
    RK_MARKER_START,
    RK_MARKER_END,
    RK_MARKER_STROKE_WIDTH,

    RK_SEL_SYM,
    RK_SEL_WIDTH,
    RK_SEL_OFFSET,
    RK_SEL_COLOR,
    RK_SEL_OPACITY,
    RK_SEL_SYM_CENTER_ANGLE,
    RK_SEL_BG_COLOR,
    RK_SEL_BG_OPACITY,
    RK_SEL_BG_WIDTH,
    RK_SEL_BG_OFFSET, 
    RK_SEL_BG_ANGLE_DELTA,
    RK_SEL_ROUNDED,

    RK_HIGHLIGHT_INTENSITY,
    RK_SHADE_INTENSITY,
    RK_BEVEL_WIDTH,
    RK_BEVEL_HIGHLIGHT_INTENSITY,
    RK_BEVEL_SHADE_INTENSITY,

    RK_SHADOW_TYPE,
    RK_SHADOW_COLOR,
    RK_SHADOW_INTENSITY,
    RK_SHADOW_LIGHT_SOURCE,
    RK_SHADOW_MAX_BLUR_LAYERS
  }

  enum knobMarkerTypes{
    RK_MT_NONE,
    RK_MT_LINE,
    RK_MT_ROUNDED,
    RK_MT_ROUNDED_FILLED,
    RK_MT_CIRCLE,
    RK_MT_CIRCLE_FILLED,
    RK_MT_CENTER_CIRCLE_TO_EDGE,
    RK_MT_CENTER_CIRCLE_TO_EDGE_FILLED
  }

  // a knob object
  class Knob : CanvasWidget {
    
    double nval = 0.5; // normalized value of param (from 0 to 1)
    double nvalc; // centered around 0 (-0.5 .. 0.5) (calculated)

    double size;
    double size_tenth, size_fortys, size_half;

    double angle_start = -135, angle_end = 135;
    double angle_width, angle_center, angle_deg, angle_rad; // calculated

    array<string> ar; // keep init params here

    // selection params
    double selBgR, selBgG, selBgB, selBgA;
    double selR, selG, selB, selA;

    double sel_w, sel_offset, sel_bg_w, sel_bg_offset, sel_rounded;
    double sel_bg_angle_delta, sel_angle_start;
    double sel_half_w;
    double rounded_w_prop;

    double s_dxp1,s_dyp1,s_dxp2,s_dyp2,s_dxp3,s_dyp3,s_dxp4,s_dyp4;

    double sbg_a_start, sbg_a_end, sbg_a_startR, sbg_a_endR;
    double sel_bg_half_w, sbg_dxp1,sbg_dyp1,sbg_dxp2,sbg_dyp2,sbg_dxp3,sbg_dyp3,sbg_dxp4,sbg_dyp4;

    double sbg_dx1, sbg_dx2, sbg_dy1, sbg_dy2;
    double sel_a_start, sel_a_end, sel_a_startR, sel_a_endR, sel_dx1, sel_dy1, sel_dx2, sel_dy2;
    double sel_min_angle_half, sel_radius;
    double round_angle_adjust, round_angle_adjustR;
    bool has_selection_bg, has_selection;
    bool sel_is_symmetric;
   
    // thumb params
    KnobThumb thumb; // thumb object
    int thumbType = 0;

    // knob constructor
    Knob(double w, double h, string params){
      // remember canvas size
      this.cw = w;
      this.ch = h;
      cw2=cw/2.0;
      ch2=ch/2.0;
      
      this.size = w; if (h>w) size = h;

      // parse params
      ar = params.split(";");

      // padding
      padl = f(ar[RK_PAD_LEFT_PX]);
      padr = f(ar[RK_PAD_RIGHT_PX]);
      padt = f(ar[RK_PAD_TOP_PX]);
      padb = f(ar[RK_PAD_BOTTOM_PX]);
      // full canvas width
      cfw = cw + padl + padr;
      cfh = ch + padt + padb;

      // position of canvas top left about center of gui
      if (ar[RK_GUI_X] == "") gui_canvas_x = -cfw/2; else gui_canvas_x = f(ar[RK_GUI_X]);
      if (ar[RK_GUI_Y] == "") gui_canvas_y = -cfh/2; else gui_canvas_y = f(ar[RK_GUI_Y]);
      
      // get canvas left top position if it is given as center position
      if (ar[RK_GUI_XY_CENTERED] == "1") {
        gui_canvas_x -= cw2;
        gui_canvas_y -= ch2;
      }

      // status("aa" + rint(gui_canvas_x)+ ", "+rint(gui_canvas_y));

      thumb.init(cw, ch, gui_canvas_x, gui_canvas_y, (ar[RK_BODY_ENABLE_3D] == "1"));

      size_tenth = size/10.0;
      size_fortys = size/40.0;
      size_half = size/2.0;

      double thumb_size = size*f(ar[RK_BODY_SIZE_RATIO]);

      // tallness and elevation depends on main size
      double thumb_tallness = thumb_size*f(ar[RK_BODY_TALLNESS_RATIO]);

      double thumb_elevation = thumb_size*f(ar[RK_BODY_ELEVATION_RATIO]);
      double thumb_round = f(ar[RK_BODY_ROUND]);
      double thumb_rotation = 0;
      int thumb_corners = parseInt(ar[RK_BODY_CORNERS]); // 4 for perfect circle :)

      thumb.setShape(thumb_corners, thumb_size, thumb_size, thumb_tallness, thumb_elevation, thumb_round, thumb_rotation);
      thumb.setBodyColor(ar[RK_BODY_COLOR], f(ar[RK_BODY_OPACITY]));

      

      // add shade and highlight
      if ((ar[RK_SHADE_INTENSITY] != "0") or (ar[RK_HIGHLIGHT_INTENSITY] != "0")) {
        thumb.setShadeAndHighlight(f(ar[RK_SHADE_INTENSITY]), f(ar[RK_HIGHLIGHT_INTENSITY]));
      }
      // add decorated bevel
      if (((ar[RK_BEVEL_SHADE_INTENSITY] != "0") or (ar[RK_BEVEL_HIGHLIGHT_INTENSITY] != "0")) and (ar[RK_BEVEL_WIDTH] != "0")) {
        thumb.setDecorBevel(f(ar[RK_BEVEL_WIDTH]), f(ar[RK_BEVEL_HIGHLIGHT_INTENSITY]), f(ar[RK_BEVEL_SHADE_INTENSITY]));
      } 
      thumb.calcShadingAndBevel();

      // set center position
      thumb.setCenterPos(cw2+padl, ch2+padt);

      float thumb_border_width = size_half * f(ar[RK_BODY_BORDER_WIDTH]);
      thumb.setBorder(ar[RK_BODY_BORDER_COLOR], f(ar[RK_BODY_BORDER_OPACITY]), thumb_border_width);

      // add marker
      thumb.setMarkerColors( ar[RK_MARKER_COLOR], f(ar[RK_MARKER_OPACITY]), "#000000", 0);
      thumb.setMarkerForKnob( parseInt(ar[RK_MARKER_TYPE]), size_tenth * f(ar[RK_MARKER_SIZE_RATIO]), (size_half - size_half * f(ar[RK_MARKER_START])), (size_half - size_half * f(ar[RK_MARKER_END])), f(ar[RK_MARKER_STROKE_WIDTH]) );

      // add shadow (shadow type, color, opacity, supported light source (-1 for all sources), max_blur_layers)
      if ((ar[RK_SHADOW_INTENSITY] != "0")) {
        string lss = ar[RK_SHADOW_LIGHT_SOURCE];
        if ((lss == "") or (lss.findFirst("all")>-1)) lss = "-1"; // all light sources if not selected
        array<string> arlss = lss.split(",");
        for(uint sn=0;sn<arlss.length;sn++) {
          if (arlss[sn].findFirst("a")>-1) arlss[sn] = "0";
          if (arlss[sn].findFirst("b")>-1) arlss[sn] = "1";
          if (arlss[sn].findFirst("c")>-1) arlss[sn] = "2";
          thumb.addShadow(parseInt(ar[RK_SHADOW_TYPE]), ar[RK_SHADOW_COLOR], f(ar[RK_SHADOW_INTENSITY]), parseInt(arlss[sn]), parseInt(ar[RK_SHADOW_MAX_BLUR_LAYERS]));
        }
      } 
      thumb.calcShadows();



      // calculate working angles
      angle_start = f(ar[RK_ANGLE_START]);
      angle_end = f(ar[RK_ANGLE_END]);

      if (angle_end < angle_start) {
        angle_start = -f(ar[RK_ANGLE_START]);
        angle_end = -f(ar[RK_ANGLE_END]);
      }

      angle_width = (angle_end-angle_start);
      angle_center = 90 - (angle_start+angle_end)/2 ; // center with shift (in normal grads)

      // prepare selection

      sel_bg_w = size_half * f(ar[RK_SEL_BG_WIDTH]);
      sel_bg_offset = size_half * f(ar[RK_SEL_BG_OFFSET]);
      sel_w = sel_bg_w * f(ar[RK_SEL_WIDTH]);
      sel_offset = size_half * f(ar[RK_SEL_OFFSET]);

      sel_rounded = f(ar[RK_SEL_ROUNDED]);
      sel_bg_angle_delta = f(ar[RK_SEL_BG_ANGLE_DELTA]);
      sel_angle_start = angle_start+270;
      sel_is_symmetric = (ar[RK_SEL_SYM] == "1");

      // if (sel_bg_w < 0) {
      //   sel_bg_w = sel_w;
      //   sel_bg_offset = sel_offset;
      // }

      sel_min_angle_half = f(ar[RK_SEL_SYM_CENTER_ANGLE]);
      // if (sel_rounded>0) sel_min_angle_half = sel_rounded*0.5;

      // calc selection radius once
      sel_radius = size_half - sel_bg_w*0.5 - sel_offset - sel_bg_offset;

      // selection bg color
      selBgR = double(parseInt(ar[RK_SEL_BG_COLOR].substr(1,2),16))/255.0;
      selBgG = double(parseInt(ar[RK_SEL_BG_COLOR].substr(3,2),16))/255.0;
      selBgB = double(parseInt(ar[RK_SEL_BG_COLOR].substr(5,2),16))/255.0;
      selBgA = f(ar[RK_SEL_BG_OPACITY]);

      // selection color
      selR = double(parseInt(ar[RK_SEL_COLOR].substr(1,2),16))/255.0;
      selG = double(parseInt(ar[RK_SEL_COLOR].substr(3,2),16))/255.0;
      selB = double(parseInt(ar[RK_SEL_COLOR].substr(5,2),16))/255.0;
      selA = f(ar[RK_SEL_OPACITY]);

      has_selection = ((selA > 0) and (sel_w > 0));
      has_selection_bg = ((selBgA > 0) and (sel_bg_w > 0));

      if (has_selection_bg) calcSelectionBg();

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

    // changes value of param
    void updateNVal(double nvalue){
      
      nval = nvalue;
      if (nval == UNKNOWN_VALUE) return; // nval is not set ()
      if (nval<0) nval=0;
      if (nval>1) nval=1;

      // calculate angle for current value
      // some functions use normal angle, others reversed
      // we get angle in "KnobMan" format with 0deg on top (on 90deg)
      nvalc = nval-0.5; // centered normalized value (-0.5 .. 0.5)
      angle_deg = angle_center - nvalc*angle_width;
      angle_rad = angle_deg*pi/180;

      // recalculate selection
      if (has_selection) calcSelection();

      // recalculate thumb and marker
      thumb.arad = angle_rad;
      thumb.adeg = angle_deg;
      thumb.calcMarker(); 
      // don't update thumb rotation if thumb is a perfect circle
      if (!thumb.is_perfect_circle) {
        thumb.setRotationDelta(nvalc*angle_width);
      }

      // status("hey: "+nval);
    }

    // on light change
    void updateRenderSettings() {
      thumb.updateRenderSettings();
    }

    // run once on init
    void calcSelectionBg(){
      
      // thank you, Pythagoras
      // calculate adjustment angle for rounded corner
      double c = (sel_w*0.5*sqrt(sel_rounded))/sel_radius;
      round_angle_adjustR = acos((2 - c*c)/2);
      round_angle_adjust = round_angle_adjustR*180/pi;

      sbg_a_start = (angle_center + 0.5*(angle_width)) + sel_bg_angle_delta - round_angle_adjust;
      sbg_a_end = (angle_center - 0.5*(angle_width)) - sel_bg_angle_delta + round_angle_adjust;

      sbg_a_startR = sbg_a_start*pi/180;
      sbg_a_endR = sbg_a_end*pi/180;

      sbg_dx1 = cos(sbg_a_startR)*(size_half - sel_bg_offset - (sel_bg_w/2)); // centers of arcs
      sbg_dy1 = -sin(sbg_a_startR)*(size_half - sel_bg_offset - (sel_bg_w/2));
      sbg_dx2 = cos(sbg_a_endR)*(size_half - sel_bg_offset - (sel_bg_w/2));
      sbg_dy2 = -sin(sbg_a_endR)*(size_half - sel_bg_offset - (sel_bg_w/2));

      // additional points for rounded corners
      sel_bg_half_w = sel_bg_w*0.5;
      float rounded_w_prop = sel_bg_half_w*(1-sel_rounded);

      sbg_dxp1 = sbg_dx1 + sin(pi2+sbg_a_startR)*rounded_w_prop;
      sbg_dyp1 = sbg_dy1 + cos(pi2+sbg_a_startR)*rounded_w_prop;

      sbg_dxp2 = sbg_dx1 + sin(-pi2+sbg_a_startR)*rounded_w_prop;
      sbg_dyp2 = sbg_dy1 + cos(-pi2+sbg_a_startR)*rounded_w_prop;

      sbg_dxp3 = sbg_dx2 + sin(pi2+sbg_a_endR)*rounded_w_prop;
      sbg_dyp3 = sbg_dy2 + cos(pi2+sbg_a_endR)*rounded_w_prop;

      sbg_dxp4 = sbg_dx2 + sin(-pi2+sbg_a_endR)*rounded_w_prop;
      sbg_dyp4 = sbg_dy2 + cos(-pi2+sbg_a_endR)*rounded_w_prop;
    }

    // run on every param change
    void calcSelection(){

      

      if (!sel_is_symmetric) {
        // if selection is normal (not symmetric)
        sel_a_start = (angle_center + 0.5*(angle_width-sel_rounded) - round_angle_adjust);
        double coeff = sel_rounded;
        if (nval > coeff) coeff=nval;
        sel_a_end = (angle_center - nvalc*angle_width + sel_rounded*0.5) + round_angle_adjust*coeff;
        // if (sel_a_end < sel_a_end_max) sel_a_end = sel_a_end_max;
      } else {
        // if selection is symmetric (around center)
        sel_a_start = (angle_center - nvalc*(angle_width-sel_rounded) - round_angle_adjust);
        sel_a_end = angle_center;
        if (sel_a_start < sel_a_end) sel_a_start = sel_a_end;
        if (nvalc > 0) {
          sel_a_start = angle_center;
          sel_a_end = (angle_center - nvalc*(angle_width-sel_rounded) + round_angle_adjust);
        }
      }
      if (sel_a_end > sel_a_start) sel_a_end = sel_a_start;
      sel_a_startR = sel_a_start*pi/180;
      sel_a_endR = sel_a_end*pi/180;

      sel_dx1 = cos(sel_a_startR)*sel_radius; // centers of arcs
      sel_dy1 = -sin(sel_a_startR)*sel_radius;
      sel_dx2 = cos(sel_a_endR)*sel_radius;
      sel_dy2 = -sin(sel_a_endR)*sel_radius;

      sel_half_w = sel_w/2;
      rounded_w_prop = sel_half_w*(1-sel_rounded);
      s_dxp1 = sel_dx1 + sin(pi2+sel_a_startR)*rounded_w_prop;
      s_dyp1 = sel_dy1 + cos(pi2+sel_a_startR)*rounded_w_prop;

      s_dxp2 = sel_dx1 + sin(-pi2+sel_a_startR)*rounded_w_prop;
      s_dyp2 = sel_dy1 + cos(-pi2+sel_a_startR)*rounded_w_prop;

      s_dxp3 = sel_dx2 + sin(pi2+sel_a_endR)*rounded_w_prop;
      s_dyp3 = sel_dy2 + cos(pi2+sel_a_endR)*rounded_w_prop;

      s_dxp4 = sel_dx2 + sin(-pi2+sel_a_endR)*rounded_w_prop;
      s_dyp4 = sel_dy2 + cos(-pi2+sel_a_endR)*rounded_w_prop;
    }

    

    // main function to re-draw the complex object
    void Draw(Kt::Graphics::Context@ ctx){

      // draw the selection
      if (has_selection or has_selection_bg) {
        DrawSelection(ctx);
      }
      
      // draw the thumb (and the marker on it)
      thumb.Draw(ctx);

      // draw debug
      // if (debug_mode != 0) DrawDebug(ctx);

    } // end of Draw of Knob

    void DrawSelection(Kt::Graphics::Context@ ctx){
      if (!has_selection and !has_selection_bg) return;
      ctx.transform.Translate(cw2+padl,ch2+padt);

      // draw selection background (if needed)
      if (has_selection_bg) {
        ctx.path.Clear();
        ctx.settings.set_lineWidth(sel_bg_w);
        ctx.source.SetRGBA (selBgR, selBgG, selBgB, selBgA);

        if (sel_rounded == 0) {
          // straight ends of selection bg
          ctx.path.Arc(0, 0, size*0.5 - sel_bg_w*0.5 - sel_bg_offset, sel_angle_start-sel_bg_angle_delta, sel_angle_start+angle_width+sel_bg_angle_delta);
          ctx.StrokePath();
        } else {
          // rounded ends of selection_bg
          ctx.path.Arc(sbg_dxp2, sbg_dyp2, sel_bg_half_w*sel_rounded, 360-180-sbg_a_start, 360-90-sbg_a_start);
          ctx.path.Arc(sbg_dxp1, sbg_dyp1, sel_bg_half_w*sel_rounded, 360-90-sbg_a_start, 360-sbg_a_start);
          ctx.path.Arc(0, 0, size_half - sel_bg_offset, 360-sbg_a_start, 360-sbg_a_end);
          ctx.path.Arc(sbg_dxp3, sbg_dyp3, sel_bg_half_w*sel_rounded, 360-sbg_a_end, 360+90-sbg_a_end);
          ctx.path.Arc(sbg_dxp4, sbg_dyp4, sel_bg_half_w*sel_rounded, 360+90-sbg_a_end, 360+180-sbg_a_end);
          ctx.path.ArcNegative(0, 0, size_half - sel_bg_offset - sel_bg_w, 360-sbg_a_end, 360-sbg_a_start);
          ctx.path.Close();
          ctx.FillPath();          
        }
      }

      ctx.SaveState();

      // now draw the selection itself
      if (sel_w > 0) {

        ctx.path.Clear();

        if (sel_rounded > 0) {
          // draw rounded corners
          ctx.path.Arc(s_dxp2, s_dyp2, sel_half_w*sel_rounded, 360-180-sel_a_start, 360-90-sel_a_start);
          ctx.path.Arc(s_dxp1, s_dyp1, sel_half_w*sel_rounded, 360-90-sel_a_start, 360-sel_a_start);
          ctx.path.Arc(0, 0, sel_radius+sel_half_w, 360-sel_a_start, 360-sel_a_end);

          ctx.path.Arc(s_dxp3, s_dyp3, sel_half_w*sel_rounded, 360-sel_a_end, 360+90-sel_a_end);
          ctx.path.Arc(s_dxp4, s_dyp4, sel_half_w*sel_rounded, 360+90-sel_a_end, 360+180-sel_a_end);
          ctx.path.ArcNegative(0, 0, sel_radius-sel_half_w, 360-sel_a_end, 360-sel_a_start);

          ctx.path.Close();
          // ctx.settings.set_lineWidth(1);
          // ctx.source.SetRGBA (0, 1, 0, 1);
          // ctx.StrokePath();

          //ctx.FillPath();    
          ctx.ClipWithPath();
        }

        // if selection has rounded corners
        // if (sel_rounded > 0) {
        //   // draw a mask for selection with rounded corners
        //   ctx.path.Arc(sel_dx1, sel_dy1, sel_bg_w*0.5, 360-180-sel_a_start, 360-sel_a_start);
        //   ctx.path.Arc(0, 0, size*0.5 - sel_bg_offset, 360-sel_a_start, 360-sel_a_end);
        //   ctx.path.Arc(sel_dx2, sel_dy2, (sel_bg_w/2), 360-sel_a_end, 360-180-sel_a_end);
        //   ctx.path.ArcNegative(0, 0, size/2 - sel_bg_offset-sel_bg_w, 360-sel_a_end, 360-sel_a_start);
        //   ctx.settings.set_lineWidth(0.01);
        //   ctx.ClipWithPath();
        //   ctx.path.Clear();
        // }


        // ctx.settings.set_lineWidth(sel_w);
        ctx.path.Clear();
        ctx.source.SetRGBA (selR, selG, selB, selA);
        

        ctx.path.Clear();
        ctx.settings.set_lineWidth(sel_w);
        // depending on selection type
        if (!sel_is_symmetric) {
          // draw normal selection (left to right)
          ctx.path.Arc(0, 0, sel_radius, sel_angle_start, sel_angle_start + nval*angle_width);
        } else {
          // symmetric selection
          if (nvalc >= 0) {
            ctx.path.Arc(0, 0, sel_radius, 360-angle_center - sel_min_angle_half, 360-angle_center + sel_min_angle_half + nvalc*(angle_width-sel_min_angle_half*2));
          } else {
            ctx.path.ArcNegative(0, 0, sel_radius, 360-angle_center + sel_min_angle_half, 360-angle_center - sel_min_angle_half + nvalc*(angle_width-sel_min_angle_half*2));
          }
        }
        ctx.StrokePath();
      }
      ctx.RestoreState();

      ctx.transform.Translate(-(cw2+padl),-(ch2+padt));
    } // end of DrawSelection

  } // end of Knob class


  // KnobThumb extending CanvasObject
  class KnobThumb : CanvasObject {

    int marker_type = RK_MT_LINE;
    double marker_size, marker_start, marker_end, marker_strokeWidth;
    double marker_size_orig, marker_start_orig, marker_end_orig, marker_strokeWidth_orig;
    double arad, adeg;
    double dxme, dyme, dxms, dyms, mdeg;

    void setMarkerForKnob(int m_type, double m_size, double m_start, double m_end, double m_strokew){
      marker_type = m_type;
      marker_size = m_size;
      if (marker_size < 2) marker_size += marker_size/2;
      marker_start = m_start;
      marker_end = m_end;
      marker_strokeWidth = m_strokew;
    }

    // calculate marker position
    void calcMarker(){
      dxme = cos(arad)*(marker_end*top_zoom);
      dyme = -sin(arad)*(marker_end*top_zoom);
      dxms = cos(arad)*(marker_start*top_zoom);
      dyms = -sin(arad)*(marker_start*top_zoom);
      mdeg = 270-adeg;
    }

    // draw the knob thumb
    void Draw(Kt::Graphics::Context@ ctx) {

      DrawGeneric(ctx);

      ctx.transform.Translate(cw2+xc_offset+top_x_offset,ch2+yc_offset+top_y_offset);

      DrawMarker(ctx);
      //ctx.source.SetRGBA(1,0,0,1);
      //ctx.WriteText("x: "+ gui_x );
      //  ctx.transform.Translate(cw2+xc_offset+top_x_offset,ch2+yc_offset+top_y_offset);
      //  ctx.transform.Rotate(-rotation); 
      //  ctx.transform.Scale(top_zoom, top_zoom);

      //  // draw thumb decor
      //  bool clipping_on = false;


      // // if shading is needed
      // if (shadeIntensity > 0) {
      //   // if (!clipping_on) { ctx.SaveState(); ctx.ClipWithPath(); clipping_on = true; }

      // }

      //  // undo "clip with path" if on
      //  if (clipping_on) { ctx.RestoreState(); clipping_on = false; }

      //  ctx.transform.Scale(1/top_zoom, 1/top_zoom);
      //  ctx.transform.Rotate(rotation); 
      ctx.transform.Translate(-(cw2+xc_offset+top_x_offset),-(ch2+yc_offset+top_y_offset));

      if (true) drawLightSources(ctx);
    }

    // draw knob marker
    void DrawMarker(Kt::Graphics::Context@ ctx){
      // draw the marker
      ctx.path.Clear();
      ctx.source.SetRGBA (markerR,markerG,markerB, markerA);
      ctx.settings.set_lineWidth(marker_size*top_zoom);

      // depending on marker type
      if (marker_type > 0)
      switch(marker_type) {
        case RK_MT_CIRCLE_FILLED: // marker - filled circle
          ctx.path.Arc(dxms, dyms, marker_size*0.75*top_zoom, 0.1, 0); // center y, x, radius, (angle a, angle b)
          ctx.FillPath();
          break;

        case RK_MT_CIRCLE: // marker - stroked circle
          ctx.path.Arc(dxms, dyms, marker_size*0.75*top_zoom, 0.1, 0); // center y, x, radius, (angle a, angle b)
          ctx.settings.set_lineWidth(marker_strokeWidth);
          ctx.StrokePath();
          break;

        case RK_MT_ROUNDED: // marker - rounded line
        case RK_MT_ROUNDED_FILLED: // marker - rounded line
          ctx.path.Arc(dxms, dyms, marker_size*0.5*top_zoom, mdeg, mdeg+180);
          ctx.path.Arc(dxme, dyme, marker_size*0.5*top_zoom, mdeg+180, mdeg+360);
          ctx.path.Close();
          if (marker_type == RK_MT_ROUNDED_FILLED) {
            ctx.FillPath();
          } else {
            ctx.settings.set_lineWidth(marker_strokeWidth);
            ctx.StrokePath();
          }
          break;

        case RK_MT_CENTER_CIRCLE_TO_EDGE:
        case RK_MT_CENTER_CIRCLE_TO_EDGE_FILLED: {// like legacy vector knobs
          float marker_center_angle_adjust = 55;
          ctx.settings.set_lineWidth(marker_strokeWidth);
          ctx.path.MoveTo(dxms, dyms);
          ctx.path.ArcNegative(0, 0, marker_size, mdeg+marker_center_angle_adjust, mdeg+180-marker_center_angle_adjust);
          ctx.path.Close();
          if (marker_type == RK_MT_CENTER_CIRCLE_TO_EDGE_FILLED) {
            ctx.FillPath();
          } else ctx.StrokePath();
        }
          break;

        default: // marker = line
          ctx.path.MoveTo(dxme, dyme);
          ctx.path.LineTo(dxms, dyms);
          ctx.StrokePath();
      }
    }

  } // end of class Knob Thumb
}



