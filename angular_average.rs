extern crate num;

fn log(x: num::complex::Complex<f64>, y: num::complex::Complex<f64>) -> num::complex::Complex<f64>
{
    (x.conj() * y).ln()
}

fn exp(x: num::complex::Complex<f64>, tau:num::complex::Complex<f64>) -> num::complex::Complex<f64>
{
    x*tau.exp()
}

fn iterate(points: &Vec<num::complex::Complex<f64>>, weights: &Vec<f64>, degree: i8 , estimate: num::complex::Complex64 )->num::complex::Complex64
{
    let mut w_sum:f64 = 0.0;
    let mut new_estimate = num::complex::Complex64::new(0.0,0.0);
    let q_minus_2:i8 = degree - 2;
    let mut cq : f64 = 0.0;
    for i in 0..points.len()
    {
        let tan_vec = log(estimate, points[i]);
        let dist : f64 = tan_vec.norm();
        if dist < 1e-6 {
            continue;
            //return points[i];
        }
        let compound_weight = weights[i] * dist.powi(q_minus_2.into());
        w_sum += compound_weight;
        new_estimate += compound_weight * tan_vec;
        cq  += weights[i] * dist.powi(degree.into());
    }
    println!("Cq {}", cq); // this is the weighted Lq norm sum cost function
    new_estimate /= w_sum;
    return exp(estimate,new_estimate);
}

fn angular_average(angle_weight_tuples: Vec<(f64,f64)>, degree: i8) -> f64{
    let mut points: Vec<num::complex::Complex<f64>> = Vec::new();
    let mut weights: Vec<f64> = Vec::new();
    for aw in angle_weight_tuples.iter() {
        points.push(num::complex::Complex::new(0.0,aw.0).exp());
        weights.push(aw.1);
    }

    //let mut product_op = |val:(&num::Complex<f64>,&f64)|->&num::Complex<f64> { &(val.1*val.0) };
    // guess it's quite complicated to perform a reduce.. or fold, as Rust cleverly remaps this age-old notion
    //let estimate2 : num::Complex<f64> = (points.iter().zip(weights.iter()).map(product_op).sum() / weights.iter().sum()).0;

    let mut w_sum:f64 = 0.0;
    let mut a_est = 0.0;
    for i in 0..points.len()
    {
        w_sum += weights[i];
        a_est += weights[i] * angle_weight_tuples[i].0;
    }
    a_est /= w_sum;
    a_est -= 1.;//1e-2;//perturb here intentionally
    let mut estimate = num::complex::Complex::new(a_est.cos(), a_est.sin());
    let mut iteration_number :i32 = 0;
    let mut old_est_angle = estimate.ln().im;
    println!("Init angle {}", to_deg(a_est));
    loop {
        iteration_number+=1;
        println!("iter: {}", iteration_number);
        let new_estimate = iterate(&points, &weights, degree, estimate);
        let new_est_angle = new_estimate.ln().im;
        println!("New estimate {} new angle estimate {}", new_estimate, to_deg(new_est_angle));
        estimate = new_estimate / new_estimate.norm();
        if (new_est_angle - old_est_angle).abs() < 1e-6{
            break;
        }
        if iteration_number > 100{
            break;
        }
        old_est_angle = new_est_angle;
    }
    return estimate.ln().im;
}

fn to_rad(x:f64)->f64{
    return x / 180.0 * std::f64::consts::PI;
}

fn to_deg(x:f64)->f64{
    return x / std::f64::consts::PI * 180.0;
}

fn main()
{
    let a = num::complex::Complex::new(1., 0.);
    let b = num::complex::Complex::new(0., 1.);
    println!("log: {}", log(b,2.0*a));
    println!("exp: {}", exp(a,log(b,a)));
    //println!("Angular average {}", to_deg(angular_average([(to_rad(0.0), 0.25), (to_rad(90.), 0.25), (to_rad(180.),0.25),(to_rad(270.0), 0.25)].to_vec(), 1)));
    //println!("Angular average {}", to_deg(angular_average([(to_rad(0.0), 0.93), (to_rad(90.), 0.33), (to_rad(180.),0.33)].to_vec(), 2)));
    println!("Angular average {}", to_deg(angular_average([(to_rad(30.0), 1./3.), (to_rad(45.), 1./3.), (to_rad(50.),1./3.)].to_vec(), 1)));
}
