// struct Parabola {
//     x_left: NotNaN<f32>,
//     x_right: NotNaN<f64>
// }

//    function AddParabola ( point u )
//    {
//       par = arc under point u;
//       if (par has its circle event, when it is removed from the beachline)
//          remove this event from the queue
//       new arcs a, b, c;
//       b.site = u;
//       a.site = c.site = par.site; // site of arc is a focus of arc
//       xl, xr  = left and right edges, which comes from point on par under u
//       xl is a normal to  (a.site, b.site);
//       xr is a normal to (b.site, c.site);
//       replace par by the sequence a, xl, b, xr, c
//       CheckCircleEvent(a);
//       CheckCircleEvent(c);
//    }
// fn add_parabola(site: Point, event_queue: &mut EventQueue) {
//     unimplemented!()
// }

//    function RemoveParabola ( Parabola p )
//    {
//       l = an arc left to p;
//       r = an arc on the right from p;
//       if (l or r have their Circle events) remove these events from the queue
//       s = the circumcenter between l.site, p.site and r.site
//       x = new edge, starts at s, normal to (l.site, r.site)
//       finish two neighbour edges xl, xr at point s
//       replace a sequence xl, p, xr by new edge x
//       CheckCircleEvent(l);
//       CheckCircleEvent(r);
//    }
// fn remove_parabola(parabola: Parabola, event_queue: &mut EventQueue) {
//     unimplemented!()
// }

//    function CheckCircleEvent(Parabola p)
//    {
//       l = arc on the left to p;
//       r = arc on the right to p;
//       xl, xr = edges by the p
//       when there is no l  OR no r  OR  l.site=r.site  RETURN
//       s = middle point, where xl and xr cross each other
//       when there is no s (edges go like\ /) RETURN
//       r = distance between s an p.site (radius of curcumcircle)
//       if s.y + r is still under the sweepline  RETURN
//       e = new circle event
//       e.parabola = p;
//       e.y = s.y + r;
//       add e into queue 
//    }
// fn check_circle_event(parabola: Parabola, event_queue: &mut EventQueue) {
//     let arc_left = arc_left_of_parabola(parabola);
//     let arc_right = arc_right_of_parabola(parabola);

//     match (arc_left, arc_right) {
//         (None, _) => { return; },
//         (_, None) => { return; },
//         (Arc{site: l_site, ..}, Arc{site: r_site, ..}) if l_site == r_site => { return; }
//     }
// }