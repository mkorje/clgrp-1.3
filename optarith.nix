{
  lib,
  stdenv,
  fetchFromGitHub,
  scons,
  pari,
  gmp,
  withPARI ? true,
  primes ? 1000,
}:

stdenv.mkDerivation {
  name = "optarith";
  src = fetchFromGitHub {
    owner = "maxwellsayles";
    repo = "liboptarith";
    rev = "e294773";
    sha256 = "sha256-l8AkEu7CZqu91zRzU3KFOC3W9E0uv6lLOXYdZTke4Oo=";
  };

  nativeBuildInputs = [ scons ];
  buildInputs = [ gmp ] ++ lib.optionals withPARI [ pari ];

  postPatch = ''
    ln -s . liboptarith

    substituteInPlace SConstruct timing/SConscript \
      --replace-fail "hasPari = os.path.exists('/usr/include/pari/pari.h') or \\" "hasPari = ${
        if withPARI then "True" else "False"
      }" \
      --replace-fail "os.path.exists('/usr/local/include/pari/pari.h')" ""

    substituteInPlace tests/SConscript timing/SConscript \
      --replace-fail "Program(" "env.Program(" \
      --replace-fail "CPPPATH='../..'" "CPPPATH='..'"

    substituteInPlace SConstruct \
      --replace-fail "env = Environment(CCFLAGS=ccflags," "env = Environment(ENV=os.environ, CCFLAGS=ccflags," \
      --replace-fail "CPPPATH='..'" "CPPPATH='.'"

    substituteInPlace tests/SConscript \
      --replace-fail "import os" "import os; env = Environment(ENV=os.environ)"

    substituteInPlace timing/SConscript \
      --replace-fail "import os.path" "import os.path; env = Environment(ENV=os.environ)"
  '';

  preBuild = ''
    $CXX -I. code_gen/gen_sqrtmodp.cc primes.c -o code_gen/gen_sqrtmodp -lgmp
    ./code_gen/gen_sqrtmodp ${toString primes}
  '';

  doCheck = true;
  checkPhase = ''
    runHook preCheck
    ./tests/test_u128
    ./tests/test_gcd
    runHook postCheck
  '';

  installPhase = ''
    runHook preInstall

    mkdir -p $out/lib
    install -m644 liboptarith.a $out/lib/
    install -m644 liboptarithxx.a $out/lib/

    mkdir -p $out/include/liboptarith/gcd
    install -m644 *.h $out/include/liboptarith/
    install -m644 gcd/*.h $out/include/liboptarith/gcd/

    mkdir -p $out/bin
    install -m755 timing/timegcd $out/bin/
    install -m755 timing/timepartial $out/bin/
    ${lib.optionalString withPARI "install -m755 timing/timepari $out/bin/"}

    runHook postInstall
  '';

  meta = with lib; {
    description = "Optimized arithmetic operations for 32, 64, and 128bit integers";
    homepage = "https://github.com/maxwellsayles/liboptarith";
    platforms = platforms.linux;
  };
}
