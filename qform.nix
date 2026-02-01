{
  lib,
  stdenv,
  fetchFromGitHub,
  scons,
  pari,
  gmp,
  optarith,
  withPARI ? true,
  doCheck ? false,
}:

stdenv.mkDerivation {
  pname = "qform";
  version = "1.0.4";
  src = fetchFromGitHub {
    owner = "mkorje";
    repo = "libqform";
    rev = "e51c843";
    sha256 = "sha256-vsIecgYCYbf4+YOfxq6pKDjeqQBznNc1fOwzP2C2Q6A=";
  };

  nativeBuildInputs = [ scons ];
  buildInputs = [
    gmp
    optarith
  ]
  ++ lib.optionals withPARI [ pari ];

  postPatch = ''
    ln -s . libqform

    substituteInPlace SConstruct \
      --replace "import os" "import os; DefaultEnvironment()['ENV'] = os.environ" \
      --replace "CPPPATH=['..']" "CPPPATH=['.']"

    substituteInPlace tests/SConscript timing/SConscript \
      --replace "source_files += ['../../libqform/libqform.a']" "source_files += ['../libqform.a']" \
      --replace "source_files += ['../../liboptarith/liboptarith.a']" "libs.append('optarith')" \
      --replace "cpppath = ['../..']" "cpppath = ['..']"

    substituteInPlace timing/SConscript \
      --replace-fail "hasPari = os.path.exists('/usr/include/pari/pari.h') or \\" "hasPari = ${
        if withPARI then "True" else "False"
      }" \
        --replace-fail "os.path.exists('/usr/local/include/pari/pari.h')" "" \
        --replace-fail "ccflags = ['-O3', '-DNDEBUG', '-Werror', '-Wall']" "ccflags = ['-O3', '-DNDEBUG', '-Werror', '-Wall', '-Wno-maybe-uninitialized']"
  ''
  + lib.optionalString withPARI ''
    sed -i 's/qfi(/Qfb0(/g' timing/time_qforms.c
  '';

  inherit doCheck;
  checkPhase = ''
    runHook preCheck
    ./tests/test_qforms --basic
    ./tests/test_qforms --pow
    runHook postCheck
  '';

  installPhase = ''
    runHook preInstall

    install -Dm644 libqform.a $out/lib/libqform.a

    mkdir -p $out/include/libqform/dbreps
    install -m644 *.h $out/include/libqform/
    install -m644 dbreps/*.h $out/include/libqform/dbreps/

    install -Dm755 timing/time_qforms $out/bin/time_qforms

    runHook postInstall
  '';

  meta = with lib; {
    description = "Ideal Arithmetic in Imaginary Quadratic Number Fields";
    homepage = "https://github.com/maxwellsayles/libqform";
    platforms = platforms.linux;
  };
}
