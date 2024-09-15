BEGIN {
	started = 0
}
{
	if (started == 1) {
		print $0
	}

	if (match($0, /geometry moments/)) {
		started = 1 - started
	}

}

END {
}
