<?xml version="1.0" encoding="utf-8" ?>
<DUMMY>
<TEMPLATE id="LMR_TPL_KNOB" 
  
  param_id="_" 
  size="64" scaling="1.0" gui_x="0" gui_y="0" gui_xy_centered="1"
  positions_count="201" pixel_range="128" mouse_area_ratio="1"
  zoom="1" offset="0.5" reverse="false" fine_ratio="10"
  disabled="false" disabled_opacity="0.2"
  visible="true" opacity="1" 

  body_color="#000000" body_opacity="1" body_size="0.9" 
  body_enable_3d="1" body_tallness="0.3" body_elevation="0" 
  border_color="#000000" border_opacity="0" border_width="0.05" 
  body_corners="4" body_round="1"

  marker_type="1" marker_size="1" marker_start="0.1" marker_end="0.9" 
  marker_color="#FFFFFF" marker_opacity="1" 
  marker_stroke_width="1" 
  angle_start="-135" angle_end="135" 
   
  selection_color="#AAAAFF" selection_opacity="0.0"
  selection_width="1" selection_offset="0" 
  selection_bg_color="#777777" selection_bg_opacity="0.0"
  selection_bg_width="0.2" selection_bg_offset="0"
  selection_bg_angle_delta="0"
  selection_rounded = "0"
  selection_sym="0" 
  selection_sym_min_center_angle="1" 
  
  shadow_intensity="1" shadow_color="#000000" shadow_light_source="-1" shadow_type="0" shadow_max_blur_layers="5"

  shade_intensity="0" highlight_intensity="0"
  bevel_width="0" bevel_highlight_intensity="0" bevel_shade_intensity="0" 

  
  
  pad_px="0" pad_top_px="0" pad_bottom_px="0" pad_left_px="0" pad_right_px="0"
  draw_on_canvas="" 
  
  color="" color1="" color2="" color3="" color4="" color5="" 
  opacity1="" opacity2="" opacity3="" opacity4="" opacity5="" 
>

  <!-- generate element_id if not given -->
  <LV id="WID" script="if (&quot;$id$&quot; == &quot;$&quot;+&quot;id$&quot;) return &quot;knob_$LMPR$_$LMID$&quot;; else return &quot;$id$&quot;;" value="" />
  <NEXT_LMID />
 
  <LV id="_body_color" value="$body_color$" /><LV id="_body_opacity" value="$body_opacity$" />
  <LV id="_marker_color" value="$marker_color$" /><LV id="_marker_opacity" value="$marker_opacity$" />
  <LV id="_selection_color" value="$selection_color$" /><LV id="_selection_opacity" value="$selection_opacity$" />
  <LV id="_selection_bg_color" value="$selection_bg_color$" /><LV id="_selection_bg_opacity" value="$selection_bg_opacity$" />
  <!-- replace with placeholders, if given -->
  <R index_list="$color$"><LV id="_body_color" value="$color$" /></R>
  <R index_list="$color1$"><LV id="_body_color" value="$color1$" /></R> <!-- color1 and color are aliases -->
  <R index_list="$color2$"><LV id="_marker_color" value="$color2$" /></R>
  <R index_list="$color3$"><LV id="_selection_color" value="$color3$" /></R>
  <R index_list="$color4$"><LV id="_selection_bg_color" value="$color4$" /></R>
  
  <R index_list="$opacity1$"><LV id="_body_opacity" value="$opacity1$" /></R>
  <R index_list="$opacity2$"><LV id="_marker_opacity" value="$opacity2$" /></R>
  <R index_list="$opacity3$"><LV id="_selection_opacity" value="$opacity3$" /></R>
  <R index_list="$opacity4$"><LV id="_selection_bg_opacity" value="$opacity4$" /></R>
  
  <LV id="w" formula="$size$" />
  <LV id="h" formula="$size$" />
  <LV id="w" formula="$w$*$scaling$" />
  <LV id="h" formula="$h$*$scaling$" />

  <LV id="knob_image_wh" formula="$size$*$body_size$*$mouse_area_ratio$*$scaling$" />

  <!-- canvas size can be bigger than element size using pad_px or pad_top_px, pad_bottom_px...-->
  <LV id="_pad_top_px" formula="$pad_top_px$+$pad_px$" />
  <LV id="_pad_bottom_px" formula="$pad_bottom_px$+$pad_px$" />
  <LV id="_pad_left_px" formula="$pad_left_px$+$pad_px$" />
  <LV id="_pad_right_px" formula="$pad_right_px$+$pad_px$" />

  <LV id="CANVAS_W" formula="$w$+$_pad_left_px$+$_pad_right_px$" />
  <LV id="CANVAS_H" formula="$h$+$_pad_top_px$+$_pad_bottom_px$" />
  <LV id="CANVAS_ID" value="cnv_$WID$" />
  <LV id="use_other_canvas" value="0" />
  <R index_list="$draw_on_canvas$" ><LV id="use_other_canvas" value="1" /></R>

  <!-- param with normalized value (from 0 to 1) -->
  <FPARAM id="$WID$_nval" min="0" max="1" formula="$WID$.zoom*(abs($reverse$-($param_id$-$param_id$.min)/($param_id$.max-$param_id$.min))) - ($WID$.zoom - 1)*$WID$.offset" exposed="true" />
  
    <LAYER_STACK internal_v_align="top" internal_h_align="left">
      <!-- <SVG_BG color="#cccccc" /> -->
      <!-- create object -->
      <SCRIPT script="LM::Knob obj_$WID$($w$, $h$, &quot;$gui_x$;$gui_y$;$gui_xy_centered$;$_pad_top_px$;$_pad_bottom_px$;$_pad_left_px$;$_pad_right_px$;$body_enable_3d$;$body_size$;$body_tallness$;$body_elevation$;$_body_color$;$_body_opacity$;$border_color$;$border_opacity$;$border_width$;$body_corners$;$body_round$;$angle_start$;$angle_end$;$marker_type$;$_marker_color$;$_marker_opacity$;$marker_size$;$marker_start$;$marker_end$;$marker_stroke_width$;$selection_sym$;$selection_width$;$selection_offset$;$_selection_color$;$_selection_opacity$;$selection_sym_min_center_angle$;$_selection_bg_color$;$_selection_bg_opacity$;$selection_bg_width$;$selection_bg_offset$;$selection_bg_angle_delta$;$selection_rounded$;$highlight_intensity$;$shade_intensity$;$bevel_width$;$bevel_highlight_intensity$;$bevel_shade_intensity$;$shadow_type$;$shadow_color$;$shadow_intensity$;$shadow_light_source$;$shadow_max_blur_layers$&quot;); 
        bool obj_$WID$_initialized = false;" />
      
      <R count="($use_other_canvas$ == 0)" index_name="__dd1">
        <CANVAS opacity="1" id="$CANVAS_ID$" width="$CANVAS_W$" height="$CANVAS_H$" render_script="if (!obj_$WID$_initialized) { obj_$WID$.setInitValue($WID$_nval); }
        obj_$WID$.Draw(Kt::Graphics::GetCurrentContext())" />
      </R>
      <R count="($use_other_canvas$ == 1)" index_name="__dd2">
        <LV id="CANVAS_ID" value="$draw_on_canvas$" />
        <ONLOAD script="obj_$CANVAS_ID$.addItem(@obj_$WID$);
          if (!obj_$WID$_initialized) { obj_$WID$.setInitValue($WID$_nval); } " />
      </R>

      <LAYER_STACK width="$w$" height="$h$" v_align="top" h_align="left" h_offset="$_pad_left_px$" v_offset="$_pad_top_px$" internal_v_align="center" internal_h_align="center" opacity="0">
        
        

        <!-- draft widget for debug -->
        <!--
        <WIDGET id="w_dbg_$WID$" width="$knob_image_wh$" height="$knob_image_wh$"><SVG_BG color="#770000" opacity="0.4" /></WIDGET>
        <PARAM_LINK from="co_debug" to="w_dbg_$WID$.display" formula="x>2" />
        -->
        
        <TEMPLATE_INNER_CONTENT />
        
        <INVISIBLE_PARAM_KNOB id="$WID$" param_id="$param_id$" opacity="$opacity$" width="$knob_image_wh$" height="$knob_image_wh$" positions_count="$positions_count$" visible="$visible$" pixel_range="$pixel_range$" images_count="1" image_orientation="vertical" cursor="$CUR_HAND$" ignore_mask_files="true" reverse="$reverse$" zoom="$zoom$" offset="$offset$" fine_ratio="$fine_ratio$" />
      

        <!-- if id is set or some of these attrs used, create control links  -->
        <!-- generate id if not given -->

        <LV id="ATTR_LINK_NEEDED" script="
          if ((&quot;$id$&quot; != &quot;$&quot;+&quot;id$&quot;) or
          (&quot;$visible$&quot; != &quot;true&quot;) or
          (&quot;$opacity$&quot; != &quot;1&quot;))
          return &quot;1&quot;; else return &quot;0&quot;;"
          />
        <R count="$ATTR_LINK_NEEDED$>0">
          <!-- "disabled" attribute -->
          <PARAM id="$id$.disabled" type="boolean" default="$disabled$" />

          <!-- link visible -->
          <FPARAM id="$CANVAS_ID$_visible" formula="
            if ($id$.disabled, $visible$, $WID$.visible)
            " />
          <PARAM_LINK from="$CANVAS_ID$_visible" to="$CANVAS_ID$.visible" />
          <!-- link display -->
          <FPARAM id="$CANVAS_ID$_display" formula="
            if ($id$.disabled, $CANVAS_ID$.display, $WID$.display)
            " />
          <PARAM_LINK from="$CANVAS_ID$_display" to="$CANVAS_ID$.display" />
          <!-- link opacity -->
          <FPARAM id="$CANVAS_ID$_opacity" formula="
            if ($id$.disabled, $WID$.opacity * $disabled_opacity$, $WID$.opacity)
            " />
          <PARAM_LINK from="$CANVAS_ID$_opacity" to="$CANVAS_ID$.opacity" />
          <!-- when disabled make our mouse element invisible -->
          <FPARAM id="$WID$_visible" formula="
            if ($id$.disabled, false, $WID$.visible)
            " />
          <PARAM_LINK from="$WID$_visible" to="$WID$.visible" />
        </R>


      </LAYER_STACK>

      <!-- redraw on param change -->
      <ACTION_TRIGGER event_id="$param_id$.value_changed;$WID$.zoom.value_changed;$WID$.offset.value_changed;$WID$.reverse.value_changed" script="obj_$WID$.updateNVal($WID$_nval); $CANVAS_ID$.Invalidate()" requires="$CANVAS_ID$.Invalidate" xasync="true" />

      <EXPOSED_OBJECTS object_ids="$WID$.*" />
      
      <!-- trigger light changes -->
      <ACTION_TRIGGER event_id="render_settings_changed.value_changed" script="obj_$WID$.updateRenderSettings(); $CANVAS_ID$.Invalidate(); " xasync="true" />

      <!-- special actions for editable mode only -->
      <!-- <IF condition="($SKIN_EDITABLE$) or ($IS_EXPORTED$ == false)"> </IF> -->

      
    </LAYER_STACK>
    <!-- <PARAM_TEXT_CONTROL param_id="$param_id$" width="100" /> -->

    

</TEMPLATE>

<DEFINE>
  
  <LMR_KNOB_BLACK_B base_type="LMR_TPL_KNOB" size="70" body_size="0.75" marker_type="3" marker_size="0.50" marker_start="0.47" marker_end="0.65" marker_color="#BBBBBB" selection_color="#5799F0" selection_opacity="0.80" selection_bg_opacity="0.40" selection_bg_width="0.23" selection_bg_offset="0.01" shadow_light_source="0" shade_intensity="1.00" highlight_intensity="0.30" bevel_width="0.1" bevel_highlight_intensity="0.70" bevel_shade_intensity="1.00" />
  <LMR_KNOB_BLACK_B_SYM base_type="LMR_KNOB_BLACK_B" selection_sym="1" />

  <LMR_KNOB_BLACK_A base_type="LMR_TPL_KNOB" size="70" body_color="#181817" body_size="0.60" selection_color="#F7BC3A" selection_opacity="1.00" selection_width="0.37" selection_bg_color="#5B5B5B" selection_bg_opacity="0.15" selection_bg_width="0.22" selection_bg_offset="0.02" selection_bg_angle_delta="4.5" marker_size="0.22" marker_start="0.47" marker_end="0.80" marker_color="#ACABA9" shade_intensity="1.00" bevel_width="0.10" bevel_highlight_intensity="0.70" bevel_shade_intensity="1.00" shadow_intensity="3.00" shadow_light_source="a" />
  <LMR_KNOB_BLACK_A_SYM base_type="LMR_KNOB_BLACK_A" selection_sym="1" />
  
  <LMR_KNOB_WHITE_A base_type="LMR_TPL_KNOB" size="70" body_color="#F9F9F9" body_size="0.63" selection_color="#7DCAFA" selection_opacity="1.00" selection_width="0.42" selection_bg_color="#EAEAEA" selection_bg_opacity="0.50" selection_bg_width="0.28" selection_bg_offset="0.03" selection_bg_angle_delta="4.5" marker_size="0.25" marker_start="0.50" marker_end="0.70" marker_color="#555555" marker_opacity="0.85" border_color="#777777" border_opacity="0.30" border_width="0.0" shade_intensity="0.30" bevel_width="0.07" bevel_shade_intensity="1.00" shadow_light_source="a" />  <LMR_KNOB_WHITE_A_SYM base_type="LMR_KNOB_WHITE_A" selection_sym="1" />
  
  <LMR_KNOB_WHITE_EVA base_type="LMR_TPL_KNOB" size="70" angle_start="-145" angle_end="145" body_color="#E8ECEE" body_size="0.70" body_enable_3d="0" body_tallness="0.20" marker_type="0" selection_color="#4AB3F4" selection_opacity="1.00" selection_bg_color="#D0D0D0" selection_bg_opacity="0.20" selection_bg_offset="0.09" selection_bg_angle_delta="35.0" selection_rounded="1.00" border_color="#FFFFFF" border_opacity="0.15" shade_intensity="0.15" />
  <LMR_KNOB_WHITE_EVA_SYM base_type="LMR_KNOB_WHITE_EVA" selection_sym="1"  />
  <LMR_KNOB_WHITE_EVA_DOT base_type="LMR_KNOB_WHITE_EVA" selection_opacity="0" marker_type="5" marker_color="#4ab3f4" marker_size="0.65" marker_opacity="1" marker_start="0.2" positions_count="7" />
  
  <LMR_KNOB_3D base_type="LMR_TPL_KNOB" size="54" body_size="0.8" body_tallness="0.4" body_color="#2c9fdf" marker_color="#EEEEEE" marker_start="0.2" marker_end="0.7" marker_size="0.95" shade_intensity="1" bevel_width="0.1" bevel_highlight_intensity="3" bevel_shade_intensity="1" />

  <LMR_KNOB_3D_BLACK base_type="LMR_KNOB_3D" body_color="#2c2c2c" size="56"  />

  <LMR_KNOB_FLAT base_type="LMR_TPL_KNOB" size="50" body_color="#F3F3F3" body_opacity="0.00" body_size="0.97" body_enable_3d="0" body_tallness="0.00" marker_type="7" marker_size="1.30" marker_start="0.28" marker_end="0.60" marker_color="#5799F0" selection_color="#5799F0" selection_opacity="1.00" selection_bg_color="#555555" selection_bg_opacity="1.00" selection_bg_width="0.28" selection_bg_offset="0.02" shadow_intensity="0.00" />
  <LMR_KNOB_FLAT_SYM base_type="LMR_KNOB_FLAT" selection_sym="1" />


  <LMR_KNOB_MARKER base_type="LMR_TPL_KNOB" size="60" body_opacity="0.00" body_size="1.00" body_enable_3d="0" marker_type="5" marker_size="0.80" marker_start="0.30" marker_color="#000000" marker_opacity="0.80" marker_stroke_width="0.90" selection_opacity="0.00" selection_bg_opacity="0.00" shadow_intensity="0.00" />


</DEFINE>


</DUMMY>