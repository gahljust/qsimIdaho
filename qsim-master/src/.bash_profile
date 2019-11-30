# .bash_profile

# Get the aliases and functions
if [ -f ~/.bashrc ]; then
	. ~/.bashrc
fi

# User specific environment and startup programs

PATH=$PATH:$HOME/bin

export PATH

if [ ! -f $HOME/.ssh/id_rsa ] ; then
	/opt/bin/pssc-config-sshkeys
fi

