#include "marschner_hair.h"

#include <iostream>

const int MAX_P = 3;
const float PI = 3.141592f;


/* MODEL CONSTANTS
 * TODO: After validation, these should be part of the hair class.

// Compute azimuthal logistic scale factor from $\beta_n$
const float beta_n = 0.3f;
const float beta_m = 0.3f;
const float alpha = 2.0f;

const float s = std::sqrt(PI/8.0f) * (0.265f * beta_n + 1.194f * beta_n*beta_n + 5.372f * std::pow(beta_n, 22.0f));

float v[4]; // Longitudinal Variance
float sin_2k_alpha[3];
float cos_2k_alpha[3];

bool hair_initialize = []() {
  v[0] = (0.726f * beta_m + 0.812f * beta_m * beta_m + 3.7f * std::pow(beta_m, 20.0f));
  v[0] = v[0] * v[0];
  v[1] = 0.25f * v[0];
  v[2] = 4.0f * v[0];
  v[3] = v[2];

  float sin_alpha = std::sin(deg2rad(alpha));
  sin_2k_alpha[0] = sin_alpha;
  cos_2k_alpha[0] = safe_sqrt(1.0f - sin_alpha * sin_alpha);

  for (int i = 1; i < 3; i++) {
    sin_2k_alpha[i] = 2.0f * cos_2k_alpha[i - 1] * sin_2k_alpha[i - 1];
    cos_2k_alpha[i] = (cos_2k_alpha[i - 1] * cos_2k_alpha[i - 1]) - (sin_2k_alpha[i - 1] * sin_2k_alpha[i - 1]);
  }

  return true;
}();
*/

inline float io(float x) {
  float val = 0;
  float x2i = 1;
  int64_t ifact = 1;
  int i4 = 1;
  // I0(x) \approx Sum_i x^(2i) / (4^i (i!)^2)
  for (int i = 0; i < 10; ++i) {
    if (i > 1) ifact *= i;
    val += x2i / (i4 * ifact * ifact);
    x2i *= x * x;
    i4 *= 4;
  }
  return val;
}

inline float log_io(float x) {
  if (x > 12)
    return x + 0.5f * (-std::log(2.0f * PI) + std::log(1.0f / x) + 1.0f / (8.0f * x));
  else
    return std::log(io(x));
}

/*
float M_p(float beta, float x) {
  float norm = (1.0f) / (beta * std::sqrt(2.0f * PI));
  float exp = std::exp((-x * x) / (2.0f * beta * beta));

  return norm * exp;
}
*/
static float Mp(float cos_theta_i, float cos_theta_o, float sin_theta_i,
  float sin_theta_o, float v) {
  float a = cos_theta_i * cos_theta_o / v;
  float b = sin_theta_i * sin_theta_o / v;
  float mp =
    (v <= .1)
    ? (std::exp(log_io(a) - b - 1 / v + 0.6931f + std::log(1 / (2 * v))))
    : (std::exp(-b) * io(a)) / (std::sinh(1 / v) * 2 * v);
  return mp;
}



glm::vec3 T(glm::vec3 sigma_a_prime, float cos_gamma) {
  // Equivalent expressions using (1 + cos(2x)) == 2cos^2(x)
  //return glm::exp(-2.0f * sigma_a * (1.0f + std::cos(2.0f * gamma)));
  
  //return glm::exp(-2.0f * sigma_a_prime * (2.0f * cos_gamma * cos_gamma));

  // PBRT Version
  return glm::exp(-2.0f * sigma_a_prime * (cos_gamma));
}

float exit_phi(float gamma_i, float gamma_t, int p) {
  return 2.0f * p * gamma_t - 2.0f * gamma_i + p * PI;
}

float dphi_dh(float eta_prime, float gamma_i, float h, int p) {
  float c = std::asin(1.0f / eta_prime);
  float num = ((6.0f * p * c) / PI - 2.0f) - 3.0f * (8.0f * p * c) / (PI * PI * PI) * gamma_i * gamma_i;
  float denom = std::sqrt(1.0f - h * h);

  return num / denom;
}

float logistic(float x, float s) {
  x = std::abs(x);
  float exp = std::exp(-x/s);
  return exp / (s * (1.0f + exp) * (1.0f + exp));
}

float logistic_cdf(float x, float s) {
  return 1.0f / (1.0f + std::exp(-x / s));
}

float trimmed_logistic(float x, float s, float a, float b) {
  return logistic(x, s) / (logistic_cdf(b, s) - logistic_cdf(a, s));
}

float sample_trimmed_logistic(float x, float s, float a, float b) {
  float k = logistic_cdf(b, s) - logistic_cdf(a, s);
  float y = -s * std::log(1 / (x * k + logistic_cdf(a, s)) - 1);
  return glm::clamp(y, a, b);
}

float N_p(int p, float phi, float s, float gamma_o, float gamma_t) {
  float delta_phi = phi - exit_phi(gamma_o, gamma_t, p);
  while (delta_phi > PI) delta_phi -= 2.0f * PI;
  while (delta_phi < -PI) delta_phi += 2.0f * PI;

  return trimmed_logistic(delta_phi, s, -PI, PI);
}

std::vector<Color> MarschnerHair::A_p(float h, float cos_theta_o, const Color& t) {
  std::vector<Color> A_p = std::vector<Color>(4);
  
  float cos_gamma_o = std::sqrt(1.0f - h * h);
  float cos_theta = cos_theta_o * cos_gamma_o;
  float f = fresnel(cos_theta, 1.0f, eta); // assumes air-hair interface, yielding 1.0f

  A_p[0] = glm::vec3(f);
  A_p[1] = (1.0f - f) * (1.0f - f) * t;
  A_p[2] = A_p[1] * f * t;

  A_p[3] = A_p[2] * f * t / (Color(1.0f) - t * f);

  return A_p;
}


Color MarschnerHair::eval(const glm::vec3& wo, const glm::vec3& wi, const Intersection& isect) {
  float h = -1.0f + isect.v * 2.0f;

  // Transform to a space where normal is in positive z-axis,
  // du is positive x-axis, dv is positive y-axis
  glm::vec3 wi_local = rotateVector(orientToZAxis(isect.n), wi);
  glm::vec3 wo_local = rotateVector(orientToZAxis(isect.n), wo);

  //std::cout << "LOCAL WO=" << wo_local.x << ", " << wo_local.y << ", " << wo_local.z << "\n";
  //std::cout << "LOCAL WI=" << wi_local.x << ", " << wi_local.y << ", " << wi_local.z << "\n";

  //return wi_local;

  float sin_theta_o = wo_local.x;
  float cos_theta_o = std::sqrt(1.0f - sin_theta_o * sin_theta_o);
  //float theta_o = std::asin(sin_theta_o);
  float phi_o = std::atan2(wo_local.z, wo_local.y);


  float sin_theta_i = wi_local.x;
  float cos_theta_i = std::sqrt(1.0f - sin_theta_i * sin_theta_i);
  //float theta_i = std::asin(sin_theta_i);
  float phi_i = std::atan2(wi_local.z, wi_local.y);

  float phi = phi_i - phi_o;

  float sin_theta_t = sin_theta_o / eta;
  float cos_theta_t = std::sqrt(1.0f - sin_theta_t * sin_theta_t);


  float eta_prime = safe_sqrt(eta * eta - sin_theta_o * sin_theta_o) / cos_theta_o;

  float sin_gamma_t = h / eta_prime;
  float cos_gamma_t = safe_sqrt(1.0f - sin_gamma_t * sin_gamma_t);
  float gamma_t = std::asin(sin_gamma_t);

  float gamma_o = std::asin(h);

  Color t = T(sigma_a / cos_theta_t, cos_gamma_t);
  std::vector<Color> Ap = A_p(h, cos_theta_o, t);
  for (const Color& c : Ap) {
    //std::cout << c.x << ", " << c.y << ", " << c.z << "\n";
  }

  glm::vec3 out = glm::vec3(0.0f);
  for (int p = 0; p < MAX_P; p++) {
    // Compute $\sin \thetao$ and $\cos \thetao$ terms accounting for scales
    float sin_theta_op, cos_theta_op;
    if (p == 0) {
      sin_theta_op = sin_theta_o * cos_2k_alpha[1] - cos_theta_o * sin_2k_alpha[1];
      cos_theta_op = cos_theta_o * cos_2k_alpha[1] + sin_theta_o * sin_2k_alpha[1];
    }

    // Handle remainder of $p$ values for hair scale tilt
    else if (p == 1) {
      sin_theta_op = sin_theta_o * cos_2k_alpha[0] + cos_theta_o * sin_2k_alpha[0];
      cos_theta_op = cos_theta_o * cos_2k_alpha[0] - sin_theta_o * sin_2k_alpha[0];
    }
    else if (p == 2) {
      sin_theta_op = sin_theta_o * cos_2k_alpha[2] + cos_theta_o * sin_2k_alpha[2];
      cos_theta_op = cos_theta_o * cos_2k_alpha[2] - sin_theta_o * sin_2k_alpha[2];
    }
    else {
      std::cout << "ILLEGAL P VALUE\n";
      exit(1);
    }

    //Color Np = Ap[p] / (2.0f * std::abs(dphi_dh(eta_prime, std::asin(h), h, p)));
    //out += M_p(betas[p], theta_h - alphas[p]) * Np;
    cos_theta_op = std::abs(cos_theta_op);
    Color Np = Ap[p] * N_p(p, phi, s, gamma_o, gamma_t);
    out += Mp(cos_theta_i, cos_theta_op, sin_theta_i, sin_theta_op, v[p]) * Np;
    //std::cout << "fsum= " << out.x << ", " << out.y << ", " << out.z << "\n";
  }

  out += Mp(cos_theta_i, cos_theta_o, sin_theta_i, sin_theta_o, v[3]) * Ap[3] / (2.0f * PI);
  
  //std::cout << "fsum= " << out.x << ", " << out.y << ", " << out.z << "\n";

  if (std::abs(wi_local.z) > 0) out /= std::abs(wi_local.z);
  
  // Todo: this shouldnt be necessary.
  if (glm::any(glm::isnan(out)) || glm::any(glm::isinf(out))) return glm::vec3(0.0f, 0.0f, 0.0f);

  //std::cout << "out_f= " << out.x << ", " << out.y << ", " << out.z << "\n";
  return out;
}

std::vector<float> MarschnerHair::get_ap_pdf(float h, float cos_theta_o) {
  float sin_theta_o = std::sqrt(1.0f - cos_theta_o * cos_theta_o);
  
  float sin_theta_t = sin_theta_o / eta;
  float cos_theta_t = std::sqrt(1.0f - sin_theta_t * sin_theta_t);

  float eta_prime = std::sqrt(eta * eta - sin_theta_o * sin_theta_o) / cos_theta_o;
  float sin_gamma_t = h / eta_prime;
  float cos_gamma_t = std::sqrt(1.0f - sin_gamma_t * sin_gamma_t);

  glm::vec3 t = T(sigma_a / cos_theta_t, cos_gamma_t);
  std::vector<Color> Ap = A_p(h, cos_theta_o, t);

  float total_luminance = 0.0f;
  for (int p = 0; p <= 3; p++) {
    total_luminance += luminance(Ap[p]);
  }

  std::vector<float> ap_pdf(4);
  for (int p = 0; p <= 3; p++) {
    ap_pdf[p] = luminance(Ap[p]) / total_luminance;
  }

  return ap_pdf;
}

std::optional<glm::vec3> MarschnerHair::sample(const Intersection& isect, const glm::vec3& wo) {
  float h = -1.0f + isect.v * 2.0f;

  glm::vec3 wo_local = rotateVector(orientToZAxis(isect.n), wo);

  float sin_theta_o = wo_local.x;
  float cos_theta_o = safe_sqrt(1.0f - sin_theta_o * sin_theta_o);
  float phi_o = std::atan2(wo_local.z, wo_local.y);

  std::vector<float> ap_pdf = get_ap_pdf(h, cos_theta_o);
  float threshold = randFloat();
  int p;
  for (p = 0; p < 3; ++p) {
    if (threshold < ap_pdf[p]) break;
    threshold -= ap_pdf[p];
  }

  // Account for hair scale tilt
  float sin_theta_op, cos_theta_op;
  if (p == 0) {
    sin_theta_op = sin_theta_o * cos_2k_alpha[1] - cos_theta_o * sin_2k_alpha[1];
    cos_theta_op = cos_theta_o * cos_2k_alpha[1] + sin_theta_o * sin_2k_alpha[1];
  }
  else if (p == 1) {
    sin_theta_op = sin_theta_o * cos_2k_alpha[0] + cos_theta_o * sin_2k_alpha[0];
    cos_theta_op = cos_theta_o * cos_2k_alpha[0] - sin_theta_o * sin_2k_alpha[0];
  }
  else if (p == 2) {
    sin_theta_op = sin_theta_o * cos_2k_alpha[2] + cos_theta_o * sin_2k_alpha[2];
    cos_theta_op = cos_theta_o * cos_2k_alpha[2] - sin_theta_o * sin_2k_alpha[2];
  }
  else {
    sin_theta_op = sin_theta_o;
    cos_theta_op = cos_theta_o;
  }

  float u = std::max(randFloat(), 1e-5f);
  float cos_theta = 1.0f + v[p] * std::log(u + (1.0f - u) * std::exp(-2.0f / v[p]));
  float sin_theta = safe_sqrt(1.0f - cos_theta * cos_theta);
  float cos_phi = std::cos(2.0f * PI * randFloat());
  float sin_theta_i = -cos_theta * sin_theta_op + sin_theta * cos_phi * cos_theta_op;
  float cos_theta_i = safe_sqrt(1.0f - sin_theta_i * sin_theta_i);


  float eta_prime = std::sqrt(eta * eta - sin_theta_o * sin_theta_o) / cos_theta_o;
  float sin_gamma_t = h / eta_prime;
  float gamma_t = std::asin(sin_gamma_t);
  float gamma_o = std::asin(h);

  float delta_phi = 0.0f;
  if (p < 3) {
    delta_phi = exit_phi(gamma_o, gamma_t, p) + sample_trimmed_logistic(randFloat(), s, -PI, PI);
  }
  else {
    delta_phi = 2.0f * PI * randFloat();
  }


  float phi_i = phi_o + delta_phi;
  glm::vec3 wi = glm::vec3(sin_theta_i,
                           cos_theta_i * std::cos(phi_i),
                           cos_theta_i * std::sin(phi_i));
   
  return  rotateVector(orientFromZAxis(isect.n), wi);
  
}

float MarschnerHair::pdf(const glm::vec3& wo, const glm::vec3& wi, const Intersection& isect) {
  //return glm::dot(isect.n, wi) / PI; // Cosine Hemisphere
  
  
  float h = -1.0f + isect.v * 2.0f;

  // Transform to a space where normal is in positive z-axis,
  // du is positive x-axis, dv is positive y-axis
  glm::vec3 wi_local = rotateVector(orientToZAxis(isect.n), wi);
  glm::vec3 wo_local = rotateVector(orientToZAxis(isect.n), wo);

  float sin_theta_o = wo_local.x;
  float cos_theta_o = safe_sqrt(1.0f - sin_theta_o * sin_theta_o);
  float phi_o = std::atan2(wo_local.z, wo_local.y);

  float sin_theta_i = wi_local.x;
  float cos_theta_i = std::sqrt(1.0f - sin_theta_i * sin_theta_i);
  float phi_i = std::atan2(wi_local.z, wi_local.y);

  float eta_prime = safe_sqrt(eta * eta - sin_theta_o * sin_theta_o) / cos_theta_o;
  float sin_gamma_t = h / eta_prime;
  float gamma_t = std::asin(sin_gamma_t);
  float gamma_o = std::asin(h);



  // Compute PDF for $A_p$ terms
  std::vector<float> ap_pdf = get_ap_pdf(h, cos_theta_o);

  // Compute PDF sum for hair scattering events
  float phi = phi_i - phi_o;
  float pdf = 0;
  for (int p = 0; p < 3; ++p) {
    float sin_theta_op, cos_theta_op;
    if (p == 0) {
      sin_theta_op = sin_theta_o * cos_2k_alpha[1] - cos_theta_o * sin_2k_alpha[1];
      cos_theta_op = cos_theta_o * cos_2k_alpha[1] + sin_theta_o * sin_2k_alpha[1];
    }
    else if (p == 1) {
      sin_theta_op = sin_theta_o * cos_2k_alpha[0] + cos_theta_o * sin_2k_alpha[0];
      cos_theta_op = cos_theta_o * cos_2k_alpha[0] - sin_theta_o * sin_2k_alpha[0];
    }
    else if (p == 2) {
      sin_theta_op = sin_theta_o * cos_2k_alpha[2] + cos_theta_o * sin_2k_alpha[2];
      cos_theta_op = cos_theta_o * cos_2k_alpha[2] - sin_theta_o * sin_2k_alpha[2];
    }
    else {
      std::cout << "ILLEGAL P VALUE\n";
      exit(1);
    }

    // Handle out-of-range $\cos \thetao$ from scale adjustment
    cos_theta_op = std::abs(cos_theta_op);
    pdf += Mp(cos_theta_i, cos_theta_op, sin_theta_i, sin_theta_op, v[p]) *
      ap_pdf[p] * N_p(p, phi, s, gamma_o, gamma_t);
  }
  pdf += Mp(cos_theta_i, cos_theta_o, sin_theta_i, sin_theta_o, v[3]) *
    ap_pdf[3] * (1.0f / (2.0f * PI));
  return pdf;
    
  
}

Color MarschnerHair::evalWithPdf(const glm::vec3& wo, const glm::vec3& wi, const Intersection& isect) {
  std::cout << "EvalWithPDF not implemented for hair...\n";
  exit(1);
}


MarschnerHair::MarschnerHair(const nlohmann::json& desc) {
  float em = desc.value("eumelanin", 1.3f);
  float pm = desc.value("pheomelanin", 0.0f);

  glm::vec3 c_e = glm::vec3(0.419, 0.697, 1.37);
  glm::vec3 c_p = glm::vec3(0.187, 0.4, 1.05);

  sigma_a = em * c_e + pm * c_p;
  eta = desc.value("eta", 1.55f);

  // Compute azimuthal logistic scale factor from $\beta_n$
  const float beta_n = desc.value("beta_n", 0.3f);
  const float beta_m = desc.value("beta_m", 0.3f);
  const float alpha = desc.value("alpha", 2.0f);

  s = std::sqrt(PI / 8.0f) * (0.265f * beta_n + 1.194f * beta_n * beta_n + 5.372f * std::pow(beta_n, 22.0f));

  v[0] = (0.726f * beta_m + 0.812f * beta_m * beta_m + 3.7f * std::pow(beta_m, 20.0f));
  v[0] = v[0] * v[0];
  v[1] = 0.25f * v[0];
  v[2] = 4.0f * v[0];
  v[3] = v[2];

  float sin_alpha = std::sin(deg2rad(alpha));
  sin_2k_alpha[0] = sin_alpha;
  cos_2k_alpha[0] = safe_sqrt(1.0f - sin_alpha * sin_alpha);

  for (int i = 1; i < 3; i++) {
    sin_2k_alpha[i] = 2.0f * cos_2k_alpha[i - 1] * sin_2k_alpha[i - 1];
    cos_2k_alpha[i] = (cos_2k_alpha[i - 1] * cos_2k_alpha[i - 1]) - (sin_2k_alpha[i - 1] * sin_2k_alpha[i - 1]);
  }
}
