with import (fetchTarball "https://github.com/NixOS/nixpkgs/archive/398a9481917df6eb4d61b8ba3237941ab32fdf67.tar.gz") {}; (pkgs.buildFHSUserEnv { name = "fhs";
	profile = ''
	      set -e
	      eval "$(micromamba shell hook -s bash)"
	      micromamba activate betaspiral
	      set +e

	    '';
    targetPkgs = pkgs: (builtins.concatLists [ [ micromamba ] [] [pkgs.which]]);
 }).env

