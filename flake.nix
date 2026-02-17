{
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-unstable";
    liboptarith = {
      url = "github:mkorje/liboptarith";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    libqform = {
      url = "github:mkorje/libqform";
      inputs.nixpkgs.follows = "nixpkgs";
      inputs.liboptarith.follows = "liboptarith";
    };
  };

  outputs =
    {
      self,
      nixpkgs,
      liboptarith,
      libqform,
    }:
    let
      system = "x86_64-linux";
      pkgs = nixpkgs.legacyPackages.${system};
    in
    {
      # packages.${system} = rec {
      #   optarith = pkgs.callPackage ./optarith.nix { };
      #   qform = pkgs.callPackage ./qform.nix { inherit optarith; };
      # };

      devShells.${system}.default = pkgs.mkShell {
        nativeBuildInputs = with pkgs; [
          autoconf
          automake
          libtool
          pkg-config
          rustPlatform.bindgenHook
        ];

        buildInputs = with pkgs; [
          # openmpi
          # gmp
          # pari
          # self.packages.${system}.optarith
          # self.packages.${system}.qform
          # gdb
          # valgrind
          # hwloc
          # python3
          # rustc
          # cargo
          # clippy
          # rustfmt
          # # clang
          rustc
          cargo
          clippy
          rustfmt
          llvmPackages_latest.libclang
          llvmPackages_latest.clang
          clang-tools
        ];

        # shellHook = ''
        #   cat << EOF
        #   make distclean
        #   automake --add-missing
        #   autoreconf
        #   ./configure CC=mpicc CPPFLAGS="-DDEBUG -DKEEP_FILES -DWITH_PARI" LIBS="-lpari -lm"
        #   ./configure CC=mpicc CPPFLAGS="-DKEEP_FILES -DWITH_PARI" LIBS="-lpari -lm"
        #   ./configure CC=mpicc CPPFLAGS="-DKEEP_FILE" LIBS="-lm"
        #   make
        #   EOF
        # '';
      };
    };
}
