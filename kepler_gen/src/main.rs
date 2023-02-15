use std::f64::consts::PI;

use peroxide::fuga::*;

#[allow(non_snake_case)]
fn main() {
    let GM = 1f64;
    let r = 4f64;
    let T = 1e+2;

    let v = Uniform(0.9, 1.1);
    let r_vec = v.sample(5).fmap(|x| x * r);
    let T_vec = r_vec.fmap(|x| (x / r).powf(3f64/2f64) * T);

    let u = Uniform(0f64, 2f64 * PI);
    let thetas = u.sample(5);
    let phis = u.sample(5);

    let mut df = DataFrame::new(vec![]);

    let mut ts: Vec<Vec<f64>> = vec![];
    let mut xs: Vec<Vec<f64>> = vec![];
    let mut ys: Vec<Vec<f64>> = vec![];
    let mut vxs: Vec<Vec<f64>> = vec![];
    let mut vys: Vec<Vec<f64>> = vec![];
    let mut Es: Vec<Vec<f64>> = vec![];
    let mut Ls: Vec<Vec<f64>> = vec![];
    let mut groups: Vec<Vec<usize>> = vec![];

    let mut group = 0usize;
    for i in 0 .. r_vec.len() {
        let r = r_vec[i];
        let T = T_vec[i];
        let theta = thetas[i];
        let phi = phis[i];

        let x = r * theta.cos();
        let y = r * theta.sin();
        let v = 0.8 * (2f64 * GM / r).sqrt();
        let vx = -v * phi.sin();
        let vy = v * phi.cos();

        let mut kepler = Kepler2D::new(x, y, vx, vy, GM);
        let dg = kepler.solve(1e-4 * T, T);
        let l = dg["x"].len();
        
        ts.push(dg["t"].to_vec());
        xs.push(dg["x"].to_vec());
        ys.push(dg["y"].to_vec());
        vxs.push(dg["vx"].to_vec());
        vys.push(dg["vy"].to_vec());
        Es.push(dg["E"].to_vec());
        Ls.push(dg["L"].to_vec());
        groups.push(vec![group; l]);
        group += 1;
    }

    // Flat vec
    let ts = ts.into_iter().flatten().collect();
    let xs = xs.into_iter().flatten().collect();
    let ys = ys.into_iter().flatten().collect();
    let vxs = vxs.into_iter().flatten().collect();
    let vys = vys.into_iter().flatten().collect();
    let Es = Es.into_iter().flatten().collect();
    let Ls = Ls.into_iter().flatten().collect();
    let groups = groups.into_iter().flatten().collect();

    df.push("t", Series::new(ts));
    df.push("x", Series::new(xs));
    df.push("y", Series::new(ys));
    df.push("vx", Series::new(vxs));
    df.push("vy", Series::new(vys));
    df.push("E", Series::new(Es));
    df.push("L", Series::new(Ls));
    df.push("g", Series::new(groups));

    df.print();

    df.write_nc("kepler.nc").unwrap();
}

#[allow(non_snake_case)]
#[derive(Debug, Clone)]
struct Kepler2D {
    x: f64,
    y: f64,
    vx: f64,
    vy: f64,
    GM: f64,
}

impl Kepler2D {
    #[allow(non_snake_case)]
    fn new(x: f64, y: f64, vx: f64, vy: f64, GM: f64) -> Self {
        Kepler2D { x, y, vx, vy, GM }
    }

    fn get_r(&self) -> f64 {
        (self.x.powi(2) + self.y.powi(2)).sqrt()
    }

    fn get_r3(&self) -> f64 {
        self.get_r().powi(3)
    }

    fn get_pos(&self) -> (f64, f64) {
        (self.x, self.y)
    }

    fn get_vel(&self) -> (f64, f64) {
        (self.vx, self.vy)
    }

    #[allow(non_snake_case)]
    fn calc_a(&self) -> (f64, f64) {
        let r3 = self.get_r3();
        let GM = self.GM;
        let x = -GM * self.x / r3;
        let y = -GM * self.y / r3;
        (x, y)
    }

    // Solve kepler problem with velocity-verlet algorithm
    #[allow(non_snake_case)]
    fn step(&mut self, dt: f64) -> (f64, f64, f64, f64, f64) {
        let (ax, ay) = self.calc_a();
        let (x, y, vx, vy) = (self.x, self.y, self.vx, self.vy);

        let vx1 = vx + 0.5 * dt * ax;
        let vy1 = vy + 0.5 * dt * ay;

        let x1 = x + dt * vx1;
        let y1 = y + dt * vy1;

        let k1 = Kepler2D::new(x1, y1, vx1, vy1, self.GM);
        let (ax1, ay1) = k1.calc_a();

        let vx2 = vx1 + 0.5 * dt * ax1;
        let vy2 = vy1 + 0.5 * dt * ay1;

        self.x = x1;
        self.y = y1;
        self.vx = vx2;
        self.vy = vy2;

        let E = 0.5 * (vx2.powi(2) + vy2.powi(2)) - self.GM / self.get_r();

        (x1, y1, vx2, vy2, E)
    }

    #[allow(non_snake_case)]
    fn solve(&mut self, dt: f64, t: f64) -> DataFrame {
        let mut t_step = 0f64;
        let mut t_vec = vec![];
        let mut x_vec = vec![];
        let mut y_vec = vec![];
        let mut vx_vec = vec![];
        let mut vy_vec = vec![];
        let mut E_vec = vec![];
        let mut L_vec = vec![];

        let (x, y) = self.get_pos();
        let (vx, vy) = self.get_vel();

        t_vec.push(t_step);
        x_vec.push(x);
        y_vec.push(y);
        vx_vec.push(vx);
        vy_vec.push(vy);
        E_vec.push(0.5 * (vx.powi(2) + vy.powi(2)) - self.GM / self.get_r());
        L_vec.push(x * vy - y * vx);


        while t_step < t {
            let (x, y, vx, vy, E) = self.step(dt);
            t_step += dt;

            t_vec.push(t_step);
            x_vec.push(x);
            y_vec.push(y);
            vx_vec.push(vx);
            vy_vec.push(vy);
            E_vec.push(E);
            L_vec.push(x * vy - y * vx);
        }

        let mut df = DataFrame::new(vec![]);
        df.push("t", Series::new(t_vec));
        df.push("x", Series::new(x_vec));
        df.push("y", Series::new(y_vec));
        df.push("vx", Series::new(vx_vec));
        df.push("vy", Series::new(vy_vec));
        df.push("E", Series::new(E_vec));
        df.push("L", Series::new(L_vec));

        df
    }
}


//fn solve_cubic(a:f64, b:f64) -> Vec<f64> {
//    let p = a.powi(2) / 3f64;
//    let q = b - 2f64 * a.powi(3) / 27f64;
//
//    let d = q.powi(2) / 4f64 - p.powi(3) / 27f64;
//
//    if d >= 0f64 {
//        let u = -q / 2f64 + d.sqrt();
//        let v = -q / 2f64 - d.sqrt();
//        vec![u.cbrt() + v.cbrt() + a / 3f64]
//    } else {
//        // Using Viete's trigonometric formula
//        let c = (-3f64 * q / (2f64 * p) * (3f64 / p).sqrt()).acos() / 3f64;
//        let v = (0 .. 3)
//            .map(|k| (c - 2f64 * PI * k as f64 / 3f64).cos() * 2f64 * (p / 3f64).sqrt() + a / 3f64).collect::<Vec<_>>();
//        v
//    }
//}
