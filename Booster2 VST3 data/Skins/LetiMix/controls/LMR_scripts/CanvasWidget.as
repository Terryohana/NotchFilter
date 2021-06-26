namespace LM {

// canvas for LM_CANVAS element, container to draw different objects
class Canvas{

  double cw, ch; // canvas width and height
  double gui_canvas_x, gui_canvas_y; // position of canvas about center of gui
  dictionary items; // collection of items that are drawn on this canvas

  Canvas(double width, double height, double gui_x = UNKNOWN_VALUE, double gui_y = UNKNOWN_VALUE) {
    cw = width;
    ch = height;

    // position of canvas top left about center of gui
    if (gui_x == UNKNOWN_VALUE) gui_canvas_x = -cw/2; else gui_canvas_x = gui_x;
    if (gui_y == UNKNOWN_VALUE) gui_canvas_y = -ch/2; else gui_canvas_y = gui_y;
  }

  // add new item (of CanvasWidget type) to this canvas
  void addItem(CanvasWidget@ obj){
    int itemssize = items.getSize();
    items.set("item"+itemssize, @obj);

    // itemssize = items.getSize();
    // status("size: " + size);
  }

  // canvas Draw (draw objects in items)
  void DrawItems(){
    if (items.isEmpty()) return;

    auto ctx = Kt::Graphics::GetCurrentContext();

    // draw debug info
    // if (co_debug > 0) {
    //   ctx.path.Clear();
    //   ctx.path.MoveTo(5,30);
    //   ctx.source.SetRGBA (.4,.4,.4, 1);
    //   ctx.WriteText("" +cw +","+ ch + " " + gui_canvas_x + ", " + gui_canvas_y);

    //   // draw draft border
    //   ctx.path.Clear();
    //   ctx.source.SetRGBA (1,.5,.5, 1);
    //   ctx.path.Rectangle(0, 0, cw, ch);
    //   ctx.settings.set_lineWidth(1);
    //   ctx.StrokePath();
    // }

    int qty = items.getSize();
    CanvasWidget@ obj;
    for(int n=0;n<qty;n++){
     if (items.get("item"+n, @obj)) {

      CanvasWidget@ ccobj = cast<CanvasWidget>(obj);

      double ccobj_gui_x = ccobj.gui_canvas_x;
      double ccobj_gui_y = ccobj.gui_canvas_y;

      if( ccobj !is null ) {
        ctx.transform.Translate(ccobj_gui_x-gui_canvas_x,ccobj_gui_y-gui_canvas_y);
        ccobj.Draw(ctx);
        ctx.transform.Translate(-1*(ccobj_gui_x-gui_canvas_x),-1*(ccobj_gui_y-gui_canvas_y));
      }
     }
    }

  }
} // end of Canvas class


// a placeholder for objects like Slider, Knob, Joystic etc
class CanvasWidget{
  int debug_mode = 0;

  float cw, ch; // canvas "working" width excluding padding
  float cfw, cfh; // canvas real width including padding (calculated)
  float padt, padb, padl, padr; // padding values
  float cw2, ch2; // half of canvas size (calculated)

  // gui offset (position of top left edge of this canvas about center of GUI)
  float gui_canvas_x, gui_canvas_y;

  // placeholder for Draw function
  void Draw(Kt::Graphics::Context@ ctx) {
    int d = 1;
  }

  // simple debug helper
  void DrawDebug(Kt::Graphics::Context@ ctx) {
    
    if (debug_mode == 0) return;

    if (debug_mode > 1) {
    ctx.path.Clear();
    ctx.path.MoveTo(5,10);
    ctx.source.SetRGBA (.4,.4,.4, 1);
    ctx.WriteText("" +gui_canvas_x +","+ gui_canvas_y);
    }

    // draw draft border
    ctx.path.Clear();
    ctx.source.SetRGBA (.5,.5,.5, 1);
    ctx.path.Rectangle(0, 0, cfw, cfh);
    ctx.settings.set_lineWidth(1);
    ctx.StrokePath();
  }

} // end of CanvasWidget class

}