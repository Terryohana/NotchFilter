<?xml version="1.0" encoding="utf-8" ?>
<DUMMY>
<!-- orientation_num can be empty or -1 (auto), 0 (vert), 1 (horiz), 2 (vert reversed,-3 (horiz reversed) -->
<TEMPLATE id="LMR_TPL_SLIDER" 
  param_id="_" visible="true" opacity="1" 
  positions_count="201" zoom="1" offset="0.5" fine_ratio="10" reverse="" 
  disabled="false" disabled_opacity="0.2"
  use_as_button="0" button_reverse="false"

  width="40" height="100" 
  scaling="1" gui_x="" gui_y="" gui_xy_centered="0" 
  margin_before="0.05" margin_after="0.05" margin_before_px="0" margin_after_px="0" 
  pad_px="0" pad_top_px="0" pad_bottom_px="0" pad_left_px="0" pad_right_px="0"
  orientation_num="" orientation="" 

  selection_color="#00DD00" selection_opacity="1" selection_bg_color="#777777" selection_bg_opacity="1" selection_width="1" selection_width_px="-1" selection_w_offset="0" selection_length="1" selection_l_offset="0" selection_round="0" selection_sym="0" selection_border_width="0" selection_border_color="" selection_border_opacity="0"
  
  thumb_color="#FFFFFF" thumb_opacity="1" thumb_type="1" 
  thumb_size="1" thumb_size_px="" thumb_size_b="1" thumb_leg_size="0" thumb_round="0"
  thumb_enable_3d="0" thumb_tallness="1" thumb_elevation="0" 
  thumb_border_color="#777777" thumb_border_opacity="1" thumb_border_width="0"
  thumb_marker_length="0" thumb_marker_thickness="1" thumb_marker_color="#FFFFFF" thumb_marker_opacity="1" thumb_marker_round="0"
  thumb_decor_type="0" thumb_decor_opacity="1" thumb_shade_intensity="0" 

  shadow_intensity="0" shadow_color="#000000" shadow_light_source="-1" shadow_type="0" shadow_max_blur_layers="10"
  
  draw_on_canvas="" debug_mode="0" extra_style=""
  use_animation="0" animation_duration_ms="80" animation_refresh_time_ms="20"

  color="" color1="" color2="" color3="" color4="" color5="" 
  opacity1="" opacity2="" opacity3="" opacity4="" opacity5="" 
  >

  <!-- generate element_id if not given -->
  <LV id="WID" script="if (&quot;$id$&quot; == &quot;$&quot;+&quot;id$&quot;) return &quot;slider_$LMPR$_$LMID$&quot;; else return &quot;$id$&quot;;" value="" />
  <NEXT_LMID />

  <LV id="w" formula="$width$*$scaling$" />
  <LV id="h" formula="$height$*$scaling$" />
  
  <LV id="_thumb_color" value="$thumb_color$" /><LV id="_thumb_opacity" value="$thumb_opacity$" />
  <LV id="_thumb_marker_color" value="$thumb_marker_color$" /><LV id="_thumb_marker_opacity" value="$thumb_marker_opacity$" />
  <LV id="_selection_color" value="$selection_color$" /><LV id="_selection_opacity" value="$selection_opacity$" />
  <LV id="_selection_bg_color" value="$selection_bg_color$" /><LV id="_selection_bg_opacity" value="$selection_bg_opacity$" />
  <LV id="_selection_border_color" value="$selection_border_color$" /><LV id="_selection_border_opacity" value="$selection_border_opacity$" />
  <!-- replace with placeholders, if given -->
  <R index_list="$color$"><LV id="_thumb_color" value="$color$" /></R>
  <R index_list="$color1$"><LV id="_thumb_color" value="$color1$" /></R> <!-- color1 and color are aliases -->
  <R index_list="$color2$"><LV id="_thumb_marker_color" value="$color2$" /></R>
  <R index_list="$color3$"><LV id="_selection_color" value="$color3$" /></R>
  <R index_list="$color4$"><LV id="_selection_bg_color" value="$color4$" /></R>
  <!-- <R index_list="$color5$"><LV id="_selection_border_color" value="$color5$" /></R> -->
  
  <R index_list="$opacity1$"><LV id="_thumb_opacity" value="$opacity1$" /></R>
  <R index_list="$opacity2$"><LV id="_thumb_marker_opacity" value="$opacity2$" /></R>
  <R index_list="$opacity3$"><LV id="_selection_opacity" value="$opacity3$" /></R>
  <R index_list="$opacity4$"><LV id="_selection_bg_opacity" value="$opacity4$" /></R>
  <!-- <R index_list="$opacity5$"><LV id="_selection_border_opacity" value="$opacity5$" /></R> -->

  <LV id="param_slider_length" formula="$h$" />
  <LV id="param_slider_orientation" value="vertical" />
  <LV id="param_slider_reverse" value="" />
  <LV id="is_horizontal" value="0" />
  <LV id="is_reversed" value="0" />

  <LV id="param_slider_width" value="$w$" />
  <LV id="param_slider_height" value="$h$" />
  <LV id="_thumb_size_px" formula="$thumb_size$*$w$" />

  <!-- check if orientation_num is set -->
  <LV id="orient" value="-1" /><R index_list="$orientation_num$"><LV id="orient" value="$orientation_num$" /></R>
  <!-- auto orientation depending on width and height -->
  <R count="($orient$ &lt; 0)*($w$>=$h$)"><LV id="orient" value="1" /></R>
  <R count="($orient$ &lt; 0)*($w$&lt;$h$)"><LV id="orient" value="0" /></R>
  <!-- if orientation parameter is used-->
  <R index_list="$orientation$"><LV id="orient_vertical" value="0" /><LV id="orient_horizontal" value="1" /><LV id="orient" value="$orient_$index$$" /></R>
  <!-- if reverse parameter is used -->
  <R index_list="$reverse$"><R count="($orient$ == 0)" index_name="_do"><LV id="orient" value="2" /></R></R>
  <R index_list="$reverse$"><R count="($orient$ == 1)" index_name="_do"><LV id="orient" value="3" /></R></R>
  <!-- horizontal -->
  <R count="($orient$ == 1)">
    <LV id="is_horizontal" value="1" />
    <LV id="param_slider_length" formula="$w$" />
    <LV id="param_slider_orientation" value="horizontal" />
    <LV id="_thumb_size_px" formula="$thumb_size$*$h$" />
  </R>
  <!-- vertical reversed -->
  <R count="($orient$ == 2)">
    <LV id="param_slider_reverse" value="true" />
    <LV id="is_reversed" value="1" />
  </R>
  <!-- horizontal reversed -->
  <R count="($orient$ == 3)">
    <LV id="is_horizontal" value="1" />
    <LV id="is_reversed" value="1" />
    <LV id="param_slider_length" formula="$w$" />
    <LV id="param_slider_orientation" value="horizontal" />
    <LV id="param_slider_reverse" value="true" />
    <LV id="_thumb_size_px" formula="$thumb_size$*$h$" />
  </R>
  
  <!-- if thumb size is set as fixed px -->
  <R index_list="$thumb_size_px$"><LV id="_thumb_size_px" formula="$thumb_size_px$" /></R>

  <!-- compute margins -->
  <LV id="m_before" formula="$margin_before$*$_thumb_size_px$ + $margin_before_px$*$scaling$" />
  <LV id="m_after" formula="$margin_after$*$_thumb_size_px$ + $margin_after_px$*$scaling$" />

  <!-- length scaling and offset -->
  <LV id="param_slider_length" formula="$param_slider_length$*$selection_length$" />
  <LV id="param_slider_pixel_range" formula="$param_slider_length$-$m_before$-$m_after$" />
  <LV id="param_slider_v_offset" value="0" />
  <LV id="param_slider_h_offset" value="0" />
  <R count="($is_horizontal$ == 0)">
    <LV id="param_slider_height" formula="$param_slider_height$*$selection_length$" />
    <LV id="param_slider_v_offset" formula="$selection_l_offset$*($h$-$param_slider_height$)/2" />
  </R>
  <R count="($is_horizontal$ == 1)">
    <LV id="param_slider_width" formula="$param_slider_width$*$selection_length$" />
    <LV id="param_slider_h_offset" formula="$selection_l_offset$*($w$-$param_slider_width$)/2" />
  </R>

  <!-- canvas size can be bigger than element size using pad_top_px, pad_bottom_px...-->
  <LV id="_pad_top_px" formula="$pad_top_px$+$pad_px$" />
  <LV id="_pad_bottom_px" formula="$pad_bottom_px$+$pad_px$" />
  <LV id="_pad_left_px" formula="$pad_left_px$+$pad_px$" />
  <LV id="_pad_right_px" formula="$pad_right_px$+$pad_px$" />

  <LV id="CANVAS_W" formula="$w$+$_pad_left_px$+$_pad_right_px$" />
  <LV id="CANVAS_H" formula="$h$+$_pad_top_px$+$_pad_bottom_px$" />
  <LV id="CANVAS_ID" value="cnv_$WID$" />
  <LV id="use_other_canvas" value="0" />
  <R index_list="$draw_on_canvas$" ><LV id="use_other_canvas" value="1" /></R>

  <!-- not used since INVISIBLE_ implemented ==>
  <LV id="image_param_slider_scaling" formula="($CANVAS_W$ > $CANVAS_H$)*$CANVAS_W$ + ($CANVAS_H$ > $CANVAS_W$)*$CANVAS_H$" />
  -->

  <!-- param with normalized value (from 0 to 1) -->
  <FPARAM id="$WID$_nval" min="0" max="1" formula="$WID$.zoom*(abs($param_slider_reverse$-($param_id$-$param_id$.min)/($param_id$.max-$param_id$.min))) - ($WID$.zoom - 1)*$WID$.offset" exposed="true" />
  
  <!-- check if use animation (mostly for switches) -->
  <R count="$use_animation$">
    <PARAM id="$WID$_nval_animation" min="0" max="1" persistent="true" />
    <PARAM_ANIMATOR id="$WID$_nval_animator" param_id="$WID$_nval_animation" refresh_time_ms="$animation_refresh_time_ms$" duration_ms="$animation_duration_ms$" persistent="true" />
    <PARAM_LINK from="$WID$_nval" to="$WID$_nval_animator" normalized="true" />
  </R>
  
  <LAYER_STACK internal_v_align="top" internal_h_align="left">
    <!-- <SVG_BG color="#cccccc" /> -->
    <!-- create slider object -->
    <SCRIPT script="LM::Slider obj_$WID$($w$, $h$, &quot;$gui_x$;$gui_y$;$gui_xy_centered$;$_pad_top_px$;$_pad_bottom_px$;$_pad_left_px$;$_pad_right_px$;$is_horizontal$;$is_reversed$;$_selection_color$;$_selection_opacity$;$_selection_bg_color$;$_selection_bg_opacity$;$selection_width$;$selection_width_px$;$selection_w_offset$;$selection_length$;$selection_l_offset$;$selection_round$;$selection_sym$;$selection_border_width$;$_selection_border_color$;$_selection_border_opacity$;$m_before$;$m_after$;$_thumb_color$;$_thumb_opacity$;$thumb_type$;$_thumb_size_px$;$thumb_size_b$;$thumb_tallness$;$thumb_elevation$;$thumb_leg_size$;$thumb_border_color$;$thumb_border_opacity$;$thumb_border_width$;$thumb_round$;$_thumb_marker_color$;$_thumb_marker_opacity$;$thumb_marker_length$;$thumb_marker_thickness$;$thumb_marker_round$;$thumb_decor_type$;$thumb_decor_opacity$;$thumb_shade_intensity$;$thumb_enable_3d$;$shadow_type$;$shadow_color$;$shadow_intensity$;$shadow_light_source$;$shadow_max_blur_layers$;$use_as_button$;$extra_style$;$debug_mode$&quot;); 
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
      <!-- <WIDGET id="w_dbg_$WID$" width="$param_slider_width$" height="$param_slider_height$" v_offset="$param_slider_v_offset$" h_offset="$param_slider_h_offset$"><SVG_BG color="#770000" opacity="0.4" /></WIDGET> -->
      
      <!-- use as regular slider -->
      <R count="($use_as_button$ == 0)">
        <INVISIBLE_PARAM_SLIDER id="$WID$" param_id="$param_id$" width="$param_slider_width$" height="$param_slider_height$" opacity="$opacity$" pixel_range="$param_slider_pixel_range$" thumb_image="" cursor="$CUR_HAND$" ignore_mask_files="true" orientation="$param_slider_orientation$" reverse="$param_slider_reverse$" v_align="center" h_align="center" v_offset="$param_slider_v_offset$" h_offset="$param_slider_h_offset$" margin_before="$m_before$" margin_after="$m_after$" positions_count="$positions_count$" visible="$visible$" zoom="$zoom$" offset="$offset$" fine_ratio="$fine_ratio$" />

      </R>

      <!-- to use as button -->
      <R count="!($use_as_button$ == 0)">
          <INVISIBLE_PARAM_BUTTON id="$WID$" param_id="$param_id$" width="$w$" height="$h$" opacity="$opacity$" image="$LM_DIR$$CONTROLS_SUBDIR$/empty_red.png" cursor="$CUR_HAND$" ignore_mask_files="true" reverse="$button_reverse$" v_align="center" h_align="center" v_offset="$param_slider_v_offset$" h_offset="$param_slider_h_offset$" positions_count="$positions_count$" visible="$visible$" zoom="$zoom$" offset="$offset$" />
      </R>


        <PARAM_LINK from="$WID$.visible" to="$CANVAS_ID$.display" id="link_$WID$_visible" />
        <PARAM_LINK from="$WID$.display" to="$CANVAS_ID$.display" />
        <PARAM_LINK from="$WID$.opacity" to="$CANVAS_ID$.opacity" />
        <!-- additional "disabled" attribute -->
        <!--
        <PARAM id="$id$.disabled" type="boolean" default="false" />
        <PARAM_LINK from="$id$.disabled" to="$CANVAS_ID$.opacity" formula="1+$disabled_opacity$-x" />
        <PARAM_LINK from="$id$.disabled" to="link_$WID$_visible.enabled" formula="1-x" />
        <PARAM_LINK from="$id$.disabled" to="$WID$.visible" formula="1-x" />
        -->

    </LAYER_STACK>

  </LAYER_STACK>
  <!-- <PARAM_TEXT_CONTROL param_id="$param_id$" width="100" /> -->
  
  <!-- redraw slider on param change -->
  <R count="$use_animation$ == 0">
    <ACTION_TRIGGER event_id="$param_id$.value_changed" script="obj_$WID$.updateNVal($WID$_nval); $CANVAS_ID$.Invalidate()" requires="$CANVAS_ID$.Invalidate" xasync="true" />
  </R>
  <R count="$use_animation$">
    <ACTION_TRIGGER event_id="$WID$_nval_animation.value_changed" script="obj_$WID$.updateNVal($WID$_nval_animation); $CANVAS_ID$.Invalidate()" requires="$CANVAS_ID$.Invalidate;$WID$_nval_animation" xasync="false" />
  </R>
  
  <EXPOSED_OBJECTS object_ids="$WID$.*;$CANVAS_ID$.*" />
  
  <!-- trigger light changes -->
  <ACTION_TRIGGER event_id="render_settings_changed.value_changed" script="obj_$WID$.updateRenderSettings(); $CANVAS_ID$.Invalidate();" xasync="true" />

</TEMPLATE>


<!-- SLIDERS -->
<DEFINE>
  <LMR_SLIDER_EMPTY_CIRCLE base_type="LMR_TPL_SLIDER" selection_bg_color="#000000" height="120" width="30" selection_bg_opacity="1" selection_opacity="0" selection_width="0.1" thumb_type="2" thumb_round="1" thumb_opacity="0" thumb_size="0.5" thumb_leg_size="0.2" thumb_border_width="2" thumb_border_color="#000000" selection_round="1" margin_before="0.9" margin_after="0.9" />
  
  <LMR_SLIDER_EMPTY_SQUARE base_type="LMR_SLIDER_EMPTY_CIRCLE" thumb_type="2" thumb_size_b="1" thumb_round="0" selection_round="0" />

  <LMR_SLIDER_CIRCLE base_type="LMR_TPL_SLIDER" height="240" width="16" selection_width="0.20" selection_bg_color="#000000" selection_bg_opacity="0.5" selection_color="#2c9fdf" selection_opacity="1" thumb_color="#2c9fdf" thumb_type="3" thumb_round="1" thumb_opacity="1" thumb_size="0.9" thumb_leg_size="0.2" thumb_border_width="0" thumb_border_color="#000000" selection_round="0.5" margin_before="0.9" margin_after="0.9" />
  
  <LMR_SLIDER_CIRCLE_H base_type="LMR_SLIDER_CIRCLE" orientation="horizontal" width="300" height="16"  />

  <LMR_SLIDER_WEB base_type="LMR_TPL_SLIDER" height="130" width="45" margin_before="0.35" margin_after="0.35" thumb_type="1" thumb_size="0.82" thumb_size_b="0.35" thumb_leg_size="0.5" thumb_elevation="0.05" thumb_tallness="0.1" thumb_round="0.3" thumb_opacity="1" thumb_border_color="#AAAAAA" thumb_border_opacity="1" thumb_border_width="0.02" selection_round="0.15" selection_width="0.75" selection_opacity="1" selection_color="#FF9933" selection_bg_color="#aaaaaa" shadow_intensity="1" scaling="1" />
  
  <LMR_SLIDER_VAL base_type="LMR_TPL_SLIDER" height="130" width="45" margin_before="0.3" margin_after="0.3" thumb_type="1" thumb_size="0.8" thumb_size_b="0.14" thumb_leg_size="0.5" thumb_color="#7efaf8" thumb_opacity="1" selection_width="0.8" selection_opacity="0.7" selection_color="#7efaf8" selection_bg_color="#42413d" scaling="1" />

  <LMR_SLIDER_BLACK_VYNIL base_type="LMR_TPL_SLIDER"  height="250" width="45" param_id="param1" margin_before="0.5" margin_after="0.5" selection_bg_color="#000000" selection_bg_opacity="1" selection_opacity="0" selection_width="0.1" selection_w_offset="-0.05" selection_round="0.2" thumb_type="1" thumb_color="#1c1c1c" thumb_opacity="1" thumb_size="0.8" thumb_size_b="0.6" thumb_leg_size="0.3" thumb_elevation="0.1" thumb_tallness="0.3" thumb_border_width="0.01" thumb_border_color="#000000" thumb_round="0.29" thumb_marker_length="0.96" thumb_marker_thickness="1" thumb_marker_color="#FFFFFF" thumb_marker_opacity="0.9" thumb_decor_type="1" thumb_decor_opacity="1"
    thumb_shade_intensity="0.3" thumb_enable_3d="1"
    shadow_intensity="1" shadow_max_blur_layers="4" />

  <LMR_SLIDER_BLACK_VYNIL_H base_type="LMR_SLIDER_BLACK_VYNIL" orientation="horizontal" width="250" height="45" />
  
  <LMR_SLIDER_TRADITIONAL base_type="LMR_TPL_SLIDER" height="350" width="40" margin_before="2" margin_after="2" selection_bg_color="#000000" selection_bg_opacity="1" selection_opacity="0" selection_width_px="4.5" selection_round="0.25" 
  thumb_type="1" thumb_decor_type="2" 
  thumb_size_px="21" thumb_size_b="2" thumb_round="0.25"  thumb_leg_size="0.3" 
  thumb_enable_3d="1" thumb_tallness="1" thumb_elevation="0.3"
  thumb_color="#000000" thumb_shade_intensity="0.25" 
  thumb_marker_length="0.96" thumb_marker_thickness="0.6" thumb_marker_color="#FFFFFF" thumb_marker_opacity="0.9" thumb_decor_opacity="0.35"
  
  shadow_intensity="1" shadow_light_source="-1" shadow_max_blur_layers="4"
  />

  <LMR_SLIDER_TRADITIONAL_H base_type="LMR_SLIDER_TRADITIONAL" orientation="horizontal" width="350" height="40" />
</DEFINE>


<!-- SWITCHES -->
<DEFINE>
  <LMR_HSWITCH_FLAT base_type="LMR_TPL_SLIDER"
  use_as_button="1" use_animation="0" positions_count="2" orientation_num="1"
  height="22" width="34" selection_width="0.75" selection_length="0.50" selection_bg_opacity="0"  selection_opacity="0" selection_round="1" 
  selection_border_width="0.13" selection_border_color="#4a4a49" selection_border_opacity="1"
  thumb_color="#4a4a49"
  thumb_type="1" thumb_round="1" thumb_opacity="1" thumb_size="0.55" thumb_leg_size="1.45" thumb_border_width="0" thumb_border_color="#000000" margin_before="0.1" margin_after="0.1" thumb_enable_3d="0" pad_px="5" 
  shadow_intensity="0"
  />
  <LMR_VSWITCH_FLAT base_type="LMR_HSWITCH_FLAT" orientation_num="0" height="34" width="22" />

  <LMR_HSWITCH_SOLID base_type="LMR_TPL_SLIDER" debug_mode="0"
  use_as_button="1" use_animation="0" positions_count="2" orientation_num="1"
  height="23" width="34" selection_width="0.75" selection_length="0.50" selection_bg_opacity="1"  selection_bg_color="#2a2a2a" selection_opacity="0" selection_round="1" 
  selection_border_width="0.0" selection_border_color="#2a2a2a" selection_border_opacity="1"
  thumb_color="#eeeeee"
  thumb_type="1" thumb_round="1" thumb_opacity="1" thumb_size="0.55" thumb_leg_size="1.45" thumb_border_width="0" thumb_border_color="#eeeeee" margin_before="0.1" margin_after="0.1" thumb_enable_3d="1" thumb_tallness="1" pad_px="5" 
  shadow_intensity="0" shadow_max_blur_layers="5"
  scaling="1"
  />
  <LMR_VSWITCH_SOLID base_type="LMR_HSWITCH_SOLID" orientation_num="0" height="34" width="23" />

  <LMR_HSWITCH_TG base_type="LMR_TPL_SLIDER" debug_mode="0"
  use_as_button="1" use_animation="0" positions_count="2" orientation_num="1"
  height="22" width="38" selection_width="0.65" selection_length="0.65" 
  selection_bg_opacity="1"  selection_bg_color="#a8a8a8" 
  selection_opacity="1" selection_color="#40a7e3" selection_round="1" 
  selection_border_width="0" selection_border_color="" selection_border_opacity="1"
  thumb_color="#eeeeee"
  thumb_type="1" thumb_round="1" thumb_opacity="1" thumb_size="0.75" thumb_leg_size="0.5" thumb_border_width="1.2" thumb_border_color="" margin_before="0.2" margin_after="0.2" thumb_enable_3d="0" thumb_tallness="1" pad_px="5" 
  shadow_intensity="0" shadow_max_blur_layers="1" extra_style="thumb_border_auto"
  scaling="1"
  />
  <LMR_VSWITCH_TG base_type="LMR_HSWITCH_TG" orientation_num="0" height="40" width="22" />

  <LMR_HSWITCH_VG base_type="LMR_TPL_SLIDER" debug_mode="0"
  use_as_button="1" use_animation="0" positions_count="2" orientation_num="1"
  height="16" width="20" selection_width="1" selection_length="0.6" 
  selection_bg_opacity="1"  selection_bg_color="#2a2a2a" 
  selection_opacity="0" selection_round="0" 
  selection_border_width="0" selection_border_color="" selection_border_opacity="1"
  thumb_color="#eeeeee"
  thumb_round="0" thumb_opacity="1" thumb_size="0.7" thumb_size_b="0.5" thumb_leg_size="2" thumb_border_width="0"  margin_before="0" margin_after="0" thumb_enable_3d="0" thumb_tallness="1" pad_px="5" 
  shadow_intensity="0" shadow_max_blur_layers="1" 
  scaling="1"
  />
  <LMR_VSWITCH_VG base_type="LMR_HSWITCH_VG" orientation_num="0" height="20" width="16" />

  <LMR_HSWITCH_3D_WHITE base_type="LMR_TPL_SLIDER" debug_mode="0"
  use_as_button="1" use_animation="0" positions_count="2" orientation_num="1"
  height="18" width="18" 
  selection_width="0.78" selection_length="0.6" 
  selection_bg_opacity="1"  selection_bg_color="#4a4a4a" 
  selection_opacity="0" selection_round="0.0" 
  selection_border_width="0" selection_border_color="#4a4a4a" selection_border_opacity="1"
  thumb_color="#eaeaea"
  thumb_round="0.2" thumb_opacity="1" thumb_size="0.7" thumb_size_b="0.8" thumb_leg_size="1.2" thumb_border_width="0.2" thumb_border_color="#aaaaaa" thumb_border_opacity="1" margin_before="0" margin_after="0" thumb_enable_3d="1" thumb_tallness="1.2" 
  thumb_decor_type="0" thumb_decor_opacity="0" thumb_shade_intensity="1"
  pad_px="10" 
  shadow_intensity="1" shadow_max_blur_layers="6" 
  scaling="1"
  />
  <LMR_VSWITCH_3D_WHITE base_type="LMR_HSWITCH_3D_WHITE" orientation_num="0" height="18" width="18" />

  <LMR_HSWITCH_3D_DARK base_type="LMR_HSWITCH_3D_WHITE" selection_width="0.8" thumb_color="#000000" thumb_decor_opacity="0.7" thumb_size_b="1" thumb_decor_type="3" shadow_intensity="1.5" scaling="1" />
  <LMR_VSWITCH_3D_DARK base_type="LMR_HSWITCH_3D_DARK" orientation_num="0" height="18" width="18" />

  <LMR_HSWITCH_3D_WROUND base_type="LMR_TPL_SLIDER" debug_mode="0"
  use_as_button="1" use_animation="0" positions_count="2" orientation_num="1"
  height="22" width="24" 
  selection_width="0.8" selection_length="0.6" 
  selection_bg_opacity="1"  selection_bg_color="#3a3a3a" 
  selection_opacity="0" selection_round="1" 
  selection_border_width="0" selection_border_color="" selection_border_opacity="1"
  thumb_color="#eaeaea"
  thumb_round="1" thumb_opacity="1" thumb_size="0.7" thumb_size_b="1" thumb_leg_size="1.2" thumb_border_width="0"  margin_before="0" margin_after="0" thumb_enable_3d="1" thumb_tallness="1.2" 
  thumb_decor_type="4" thumb_decor_opacity="1" thumb_shade_intensity="1"
  pad_px="10" 
  shadow_intensity="1" shadow_max_blur_layers="6" 
  scaling="1"
  />
  <LMR_VSWITCH_3D_WROUND base_type="LMR_HSWITCH_3D_WROUND" orientation_num="0" height="24" width="22" />
</DEFINE>

</DUMMY>